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
    kseq_t *seq = kseq_init(fp);
    FILE *out = tmpfile();
    EXPECT_TRUE(fmsi2::run_queries<hashes>(index, seq, out, streaming, 32, threads, block_bases, stats));
    kseq_destroy(seq);
    gzclose(fp);
    std::string result(ftell(out), '\0');
    rewind(out);
    EXPECT_EQ(fread(&result[0], 1, result.size(), out), result.size());
    fclose(out);
    return result;
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
