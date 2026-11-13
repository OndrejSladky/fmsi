#pragma once

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

#include "../src/query_pipeline.h"
#include "fms_index2_test.h"
#include "gtest/gtest.h"

namespace {

/// The output of the pipeline for a FASTA file.
template <bool hashes>
std::string pipeline_output(const fmsi2::index &index, std::string path, bool streaming, unsigned threads,
                            size_t block_bases, fmsi2::query_stats &stats) {
    gzFile fp = OpenFile(path);
    FILE *out = tmpfile();
    EXPECT_TRUE(fmsi2::run_queries<hashes>(index, fp, out, streaming, 32, threads, block_bases, stats));
    gzclose(fp);
    std::string result(ftell(out), '\0');
    rewind(out);
    EXPECT_EQ(fread(&result[0], 1, result.size(), out), result.size());
    fclose(out);
    return result;
}

typedef std::vector<std::pair<std::string, std::string>> named_sequences;

/// Names and sequences read from a whole file by kseq, as by the stable commands.
named_sequences kseq_records(std::string path) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);
    named_sequences result;
    while (kseq_read(seq) >= 0) {
        result.push_back({std::string(seq->name.s, seq->name.l), std::string(seq->seq.s, seq->seq.l)});
    }
    kseq_destroy(seq);
    gzclose(fp);
    return result;
}

/// The same from the blocks of the reader of the pipeline.
named_sequences block_records(std::string path, size_t block_bytes, size_t read_step) {
    gzFile fp = OpenFile(path);
    named_sequences result;
    {
        fmsi2::block_reader reader(fp, block_bytes, read_step);
        fmsi2::query_block b;
        while (reader.read(b)) {
            if (!b.parsed) fmsi2::parse_block(b);
            for (const auto &r : b.records) {
                result.push_back({b.names.substr(r.name, r.name_length), b.sequences.substr(r.offset, r.length)});
            }
        }
    }
    gzclose(fp);
    return result;
}

std::string random_line(std::mt19937_64 &rng, const std::string &alphabet, size_t max_length) {
    std::string line(rng() % (max_length + 1), ' ');
    for (auto &c : line) c = alphabet[rng() % alphabet.size()];
    return line;
}

std::string random_fasta(std::mt19937_64 &rng) {
    std::string text;
    for (int i = rng() % 6; i > 0; --i) {
        text += ">" + random_line(rng, "abc123+>@_", 6);
        if (rng() % 2) text += " " + random_line(rng, "comment +>@", 8);
        text += rng() % 5 ? "\n" : "\r\n";
        for (int lines = rng() % 4; lines > 0; --lines) {
            text += random_line(rng, "ACGTNacgt", 30) + (rng() % 5 ? "\n" : "\r\n");
            if (rng() % 6 == 0) text += "\n";
        }
    }
    return text;
}

std::string random_fastq(std::mt19937_64 &rng) {
    std::string text;
    for (int i = rng() % 5; i > 0; --i) {
        std::string sequence = random_line(rng, "ACGTN", 30);
        std::string quality = random_line(rng, "!I5>@+#", 0);
        while (quality.size() < sequence.size()) quality += "!I5>@+#"[rng() % 7];
        text += "@" + random_line(rng, "abc123", 6) + "\n" + sequence + "\n+\n" + quality + "\n";
    }
    return text;
}

/// The reader of the pipeline cuts FASTA files into blocks without parsing them; the sequences must be
/// those that kseq reads from the whole file, for well-formed and for malformed files alike.
TEST(QUERY_PIPELINE, BLOCKS_MATCH_KSEQ) {
    std::mt19937_64 rng(12);
    std::string path = (std::filesystem::temp_directory_path() / "fmsi2_block_reader_test.txt").string();
    for (int round = 0; round < 1500; ++round) {
        std::string text;
        switch (round % 6) {
        case 0:
            text = random_fasta(rng);
            break;
        case 1:
            text = random_fastq(rng);
            break;
        case 2:
            text = random_line(rng, "junk >@+\n", 12) + random_fasta(rng) + random_fastq(rng) + random_fasta(rng);
            break;
        case 3:
            text = random_fasta(rng) + random_fasta(rng) + random_fasta(rng);
            break;
        case 4:
            // Anything, with many line starts.
            text = random_line(rng, "ACGTN>@+\n\n\n\r x", 300);
            break;
        default:
            text = random_fastq(rng) + random_fasta(rng);
            break;
        }
        // Sometimes the file ends in the middle of a record.
        if (round % 4 == 3 && !text.empty()) text.resize(rng() % text.size());
        {
            std::ofstream file(path, std::ios::binary);
            file << text;
        }
        named_sequences expected = kseq_records(path);
        // The file is read a few bytes at a time, so that records span several reads and blocks.
        for (size_t block_bytes : {(size_t)1, (size_t)3, (size_t)17, (size_t)64, (size_t)1 << 20}) {
            for (size_t read_step : {(size_t)1, (size_t)5, (size_t)1 << 20}) {
                ASSERT_EQ(block_records(path, block_bytes, read_step), expected)
                    << "round " << round << " block " << block_bytes << " step " << read_step << " text:\n" << text;
            }
        }
    }
    std::filesystem::remove(path);
}

/// The lines must follow the order of the input and agree with the definition, whatever the number
/// of threads and the size of the blocks.
TEST(QUERY_PIPELINE, OUTPUT_INDEPENDENT_OF_THREADS_AND_BLOCKS) {
    std::mt19937_64 rng(11);
    for (int round = 0; round < 6; ++round) {
        size_t k = std::vector<size_t>{5, 13, 31}[round % 3];
        std::string ms = random_mask(rng, random_superstring(rng, 3000 + rng() % 3000), k, true);
        auto represented = represented_kmers(ms, k);
        mbwt m = compute_mbwt<uint64_t>(ms, k, round % 2);
        fmsi2::index index;
        fmsi2::build(index, m, k, (fmsi2::mask_kind)(round % 3));
        fmsi2::build_prefix_table(index, 4);

        // Reads of various lengths, some shorter than k, empty or with invalid characters, and a long one.
        std::string superstring = ms, path = (std::filesystem::temp_directory_path() / "fmsi2_query_pipeline_test.fa").string();
        for (auto &c : superstring) c = toupper(c);
        std::string expected;
        {
            std::ofstream fasta(path);
            for (int i = 0; i < 400; ++i) {
                std::string read;
                if (i == 100) {
                    read = superstring + reverse_complement(superstring);
                } else if (i % 50 != 7) {
                    read = superstring.substr(rng() % superstring.size(), rng() % (6 * k));
                    if (rng() % 2) read = reverse_complement(read);
                    for (auto &c : read) {
                        if (rng() % 30 == 0) c = "ACGTN"[rng() % 5];
                    }
                }
                fasta << ">read" << i << " comment\n" << read << "\n";
                expected += "read" + std::to_string(i) + "\t";
                for (size_t j = 0; j + k <= read.size(); ++j) {
                    std::string kmer = read.substr(j, k);
                    expected += kmer.find('N') == std::string::npos && represented.count(canonical(kmer)) ? '1' : '0';
                }
                expected += "\n";
            }
        }

        fmsi2::query_stats reference_stats;
        std::string hashes = pipeline_output<true>(index, path, false, 1, 1 << 20, reference_stats);
        for (bool streaming : {false, true}) {
            for (unsigned threads : {1, 2, 5}) {
                for (size_t block_bases : {(size_t)1, (size_t)200, (size_t)1 << 20}) {
                    fmsi2::query_stats stats, hash_stats;
                    ASSERT_EQ(pipeline_output<false>(index, path, streaming, threads, block_bases, stats), expected)
                        << "round " << round << " streaming " << streaming << " threads " << threads << " block " << block_bases;
                    ASSERT_EQ(pipeline_output<true>(index, path, streaming, threads, block_bases, hash_stats), hashes)
                        << "round " << round << " streaming " << streaming << " threads " << threads << " block " << block_bases;
                    ASSERT_EQ(stats.lookups.kmers, reference_stats.lookups.kmers);
                    ASSERT_EQ(stats.lookups.found, reference_stats.lookups.found);
                }
            }
        }
        std::filesystem::remove(path);
    }
}

} // namespace
