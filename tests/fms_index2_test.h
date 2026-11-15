#pragma once

#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "../src/doubled.h"
#include "../src/fms_index2.h"
#include "../src/parser.h"
#include "gtest/gtest.h"

namespace {

std::string random_dna(std::mt19937_64 &rng, size_t length) {
    std::string result(length, 'A');
    for (auto &c : result) c = "ACGT"[rng() % 4];
    return result;
}

std::string reverse_complement(const std::string &s) {
    char *rc = ReverseComplementString(s.data(), s.size());
    std::string result(rc, s.size());
    delete[] rc;
    return result;
}

std::string canonical(const std::string &kmer) { return std::min(kmer, reverse_complement(kmer)); }

/// A random superstring with repeats on both strands, so that k-mers have several occurrences.
std::string random_superstring(std::mt19937_64 &rng, size_t length) {
    std::string s = random_dna(rng, std::min(length, (size_t)20));
    while (s.size() < length) {
        size_t piece = 1 + rng() % 40;
        switch (rng() % 3) {
        case 0:
            s += random_dna(rng, piece);
            break;
        case 1: {
            size_t from = rng() % s.size();
            s += s.substr(from, piece);
            break;
        }
        default: {
            size_t from = rng() % s.size();
            s += reverse_complement(s.substr(from, piece));
            break;
        }
        }
    }
    s.resize(length);
    return s;
}

/// Set the mask (letter case) of the superstring. With [single], each represented k-mer has exactly
/// one set occurrence (as in masks with the minimum number of ones); otherwise the bits are random.
std::string random_mask(std::mt19937_64 &rng, std::string s, size_t k, bool single) {
    for (auto &c : s) c = tolower(c);
    std::set<std::string> used;
    for (size_t i = 0; i + k <= s.size(); ++i) {
        if (rng() % 3 == 0) continue;
        std::string kmer = s.substr(i, k);
        for (auto &c : kmer) c = toupper(c);
        if (single && !used.insert(canonical(kmer)).second) continue;
        s[i] = toupper(s[i]);
    }
    return s;
}

/// The represented canonical k-mers: those with a set occurrence on either strand.
std::set<std::string> represented_kmers(const std::string &ms, size_t k) {
    std::set<std::string> result;
    for (size_t i = 0; i + k <= ms.size(); ++i) {
        if (!is_upper(ms[i])) continue;
        std::string kmer = ms.substr(i, k);
        for (auto &c : kmer) c = toupper(c);
        result.insert(canonical(kmer));
    }
    return result;
}

struct queries {
    std::string sequences;
    std::vector<uint32_t> offsets;
    std::vector<fmsi2::segment> segments;
};

queries make_queries(std::mt19937_64 &rng, const std::string &ms, size_t k) {
    std::string superstring = ms;
    for (auto &c : superstring) c = toupper(c);
    std::vector<std::string> sequences = {superstring, reverse_complement(superstring),
                                          random_dna(rng, superstring.size())};
    // Substrings in both orientations with substitutions.
    for (int i = 0; i < 30; ++i) {
        size_t from = rng() % superstring.size();
        std::string read = superstring.substr(from, k + rng() % (3 * k));
        if (rng() % 2) read = reverse_complement(read);
        for (auto &c : read) {
            if (rng() % 25 == 0) c = "ACGT"[rng() % 4];
        }
        sequences.push_back(read);
    }
    queries result;
    for (auto &sequence : sequences) {
        if (sequence.size() < k) continue;
        result.segments.push_back({(uint32_t)result.sequences.size(), (uint32_t)sequence.size()});
        for (size_t i = 0; i + k <= sequence.size(); ++i) result.offsets.push_back(result.sequences.size() + i);
        result.sequences += sequence;
    }
    return result;
}

template <size_t BYTES>
void test_quad_rank() {
    std::mt19937_64 rng(1);
    const size_t block = fmsi2::quad_rank<BYTES>::BLOCK_SYMBOLS, super = block * fmsi2::quad_rank<BYTES>::SUPER_BLOCKS;
    for (size_t size : {(size_t)0, (size_t)1, block - 1, block, block + 1, 9 * block + 5, super, super + 1, 2 * super + 77}) {
        std::vector<uint8_t> symbols(size);
        fmsi2::quad_rank<BYTES> ranks(size);
        for (size_t i = 0; i < size; ++i) {
            symbols[i] = rng() % 4;
            ranks.set(i, symbols[i]);
        }
        ranks.finalize();
        size_t counts[4] = {0, 0, 0, 0};
        for (size_t i = 0; i <= size; ++i) {
            for (int c = 0; c < 4; ++c) ASSERT_EQ(ranks.rank(i, c), counts[c]) << "size " << size << " position " << i;
            for (size_t distance : {0, 1, 2, 3, 5, 17, 27, 28, 29, 30, 250, 1100}) {
                if (i + distance > size) continue;
                size_t rank_i, rank_j;
                uint8_t c = rng() % 4;
                ranks.rank_pair(i, i + distance, c, rank_i, rank_j);
                ASSERT_EQ(rank_i, ranks.rank(i, c));
                ASSERT_EQ(rank_j, ranks.rank(i + distance, c));
                if (distance == 0) continue;
                ranks.rank_pair_near(i, i + distance, c, rank_i, rank_j);
                ASSERT_EQ(rank_i, ranks.rank(i, c));
                ASSERT_EQ(rank_j, ranks.rank(i + distance, c)) << "size " << size << " position " << i << " distance "
                                                                << distance;
                uint8_t symbol;
                ASSERT_EQ(ranks.rank_symbol(i, c, symbol), rank_i);
                if (i < size) ASSERT_EQ(symbol, symbols[i]);
            }
            if (i < size) {
                ASSERT_EQ(ranks[i], symbols[i]);
                counts[symbols[i]]++;
            }
        }
    }
}

TEST(RANK2, QUAD_RANK) {
    test_quad_rank<64>();
    test_quad_rank<128>();
    test_quad_rank<256>();
}

template <size_t BYTES>
void test_quad_select() {
    std::mt19937_64 rng(11);
    const size_t block = fmsi2::quad_rank<BYTES>::BLOCK_SYMBOLS, super = block * fmsi2::quad_rank<BYTES>::SUPER_BLOCKS;
    // Symbols drawn evenly, with long runs, and with one symbol nearly absent: the guesses of select are
    // good in the first case and bad in the others.
    for (int kind = 0; kind < 3; ++kind) {
        for (size_t size : {(size_t)1, block - 1, block, block + 1, 9 * block + 5, super, super + 1, 2 * super + 77}) {
            std::vector<uint8_t> symbols(size);
            fmsi2::quad_rank<BYTES> ranks(size);
            uint8_t current = 0;
            for (size_t i = 0; i < size; ++i) {
                if (kind == 0) {
                    current = rng() % 4;
                } else if (kind == 1) {
                    if (rng() % 700 == 0) current = rng() % 4;
                } else {
                    current = rng() % 5000 == 0 ? 3 : rng() % 3;
                }
                symbols[i] = current;
                ranks.set(i, symbols[i]);
            }
            ranks.finalize();
            // Without the select hints and with them.
            for (int hints = 0; hints < 2; ++hints) {
                if (hints) ranks.build_select_hints();
                ASSERT_EQ(ranks.has_select_hints(), (bool)hints);
                size_t counts[4] = {0, 0, 0, 0};
                for (size_t i = 0; i < size; ++i) {
                    uint8_t c = symbols[i];
                    ASSERT_EQ(ranks.select(counts[c], c), i)
                        << "kind " << kind << " size " << size << " symbol " << (int)c << " hints " << hints;
                    counts[c]++;
                }
            }
        }
    }
}

TEST(RANK2, QUAD_SELECT) {
    test_quad_select<64>();
    test_quad_select<128>();
    test_quad_select<256>();
}

TEST(RANK2, BIT_SELECT) {
    std::mt19937_64 rng(12);
    for (int kind = 0; kind < 3; ++kind) {
        for (size_t size : {1, 495, 496, 497, 5000, 63488, 63489, 140000, 400000}) {
            fmsi2::bit_rank ranks(size);
            std::vector<size_t> ones;
            bool current = true;
            for (size_t i = 0; i < size; ++i) {
                if (kind == 0) {
                    current = rng() % 3 == 0;
                } else if (kind == 1) {
                    if (rng() % 3000 == 0) current = !current;
                } else {
                    current = rng() % 20000 == 0;
                }
                if (current) {
                    ranks.set(i);
                    ones.push_back(i);
                }
            }
            ranks.finalize();
            for (size_t r = 0; r < ones.size(); ++r) {
                ASSERT_EQ(ranks.select(r), ones[r]) << "kind " << kind << " size " << size << " one " << r;
            }
        }
    }
}

TEST(RANK2, BIT_RANK) {
    std::mt19937_64 rng(2);
    for (size_t size : {0, 1, 495, 496, 497, 5000, 63488, 63489, 140000}) {
        std::vector<bool> bits(size);
        fmsi2::bit_rank ranks(size);
        for (size_t i = 0; i < size; ++i) {
            bits[i] = rng() % 3 == 0;
            if (bits[i]) ranks.set(i);
        }
        ranks.finalize();
        size_t ones = 0;
        for (size_t i = 0; i <= size; ++i) {
            ASSERT_EQ(ranks.rank(i), ones) << "size " << size << " position " << i;
            if (i < size) {
                ASSERT_EQ(ranks[i], bits[i]);
                ones += bits[i];
            }
        }
    }
}

TEST(RANK2, TUNSTALL_BITS) {
    std::mt19937_64 rng(7);
    for (double density : {0.0, 0.003, 0.228, 0.5, 0.772, 0.97, 1.0}) {
        for (size_t size : {0, 1, 2, 100, 594, 5000, 150000}) {
            std::vector<bool> bits(size);
            for (size_t i = 0; i < size; ++i) bits[i] = (rng() >> 11) * (1.0 / (1ULL << 53)) < density;
            fmsi2::tunstall_bits ranks;
            ranks.build(size, [&](size_t i) { return (bool)bits[i]; });
            size_t ones = 0;
            for (size_t i = 0; i <= size; ++i) {
                ASSERT_EQ(ranks.rank(i), ones) << "density " << density << " size " << size << " position " << i;
                if (i < size) {
                    ASSERT_EQ(ranks[i], bits[i]) << "density " << density << " size " << size << " position " << i;
                    ones += bits[i];
                }
            }
        }
    }
    // Runs, which the code is not made for, must still be decoded right.
    size_t size = 100000;
    std::vector<bool> bits(size);
    for (size_t i = 0; i < size; ++i) bits[i] = (i / 3000) % 2;
    fmsi2::tunstall_bits ranks;
    ranks.build(size, [&](size_t i) { return (bool)bits[i]; });
    size_t ones = 0;
    for (size_t i = 0; i < size; ++i) {
        ASSERT_EQ(ranks.rank(i), ones);
        ASSERT_EQ(ranks[i], bits[i]);
        ones += bits[i];
    }
}

TEST(FMS_INDEX2, MATCHES_SUFFIX_ARRAY_CONSTRUCTIONS) {
    std::mt19937_64 rng(3);
    for (int round = 0; round < 20; ++round) {
        size_t k = 2 + rng() % 30;
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 3000), k, false);
        mbwt libsais = compute_mbwt<uint64_t>(ms, k, true, false);
        mbwt qsufsort = compute_mbwt<uint64_t>(ms, k, true, true);
        EXPECT_EQ(libsais.bwt, qsufsort.bwt);
        EXPECT_EQ(libsais.sa_transformed_mask, qsufsort.sa_transformed_mask);
        EXPECT_EQ(libsais.klcp, qsufsort.klcp);
        EXPECT_EQ(libsais.dollar_position, qsufsort.dollar_position);
    }
}

/// The index must not depend on the number of threads that construct it.
TEST(FMS_INDEX2, CONSTRUCTION_INDEPENDENT_OF_THREADS) {
    std::mt19937_64 rng(6);
    // Small chunks, so that the loops over the short test strings are split among the threads.
    size_t usual_chunk = parallel_chunk_items;
    parallel_chunk_items = 64;
    auto file_bytes = [](const std::string &path) {
        std::ifstream in(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    };
    std::string path = (std::filesystem::temp_directory_path() / "fmsi2_construction_test").string();
    for (int round = 0; round < 12; ++round) {
        size_t k = std::vector<size_t>{2, 5, 13, 31}[round % 4];
        std::string ms = random_mask(rng, random_superstring(rng, k + 5000 + rng() % 60000), k, round % 2);
        mbwt single = compute_mbwt<uint64_t>(ms, k, true, false, 1);
        mbwt threaded = compute_mbwt<uint64_t>(ms, k, true, false, 3);
        ASSERT_EQ(single.bwt, threaded.bwt);
        ASSERT_EQ(single.sa_transformed_mask, threaded.sa_transformed_mask);
        ASSERT_EQ(single.klcp, threaded.klcp);
        ASSERT_EQ(single.dollar_position, threaded.dollar_position);

        auto mask_type = (fmsi2::mask_kind)(round % 3);
        fmsi2::index a, b;
        fmsi2::build(a, single, k, mask_type, 64 << (round % 3), 1);
        fmsi2::build(b, threaded, k, mask_type, 64 << (round % 3), 5);
        fmsi2::save(a, path + ".1");
        fmsi2::save(b, path + ".5");
        ASSERT_EQ(file_bytes(path + ".1"), file_bytes(path + ".5")) << "round " << round;
    }
    std::filesystem::remove(path + ".1");
    std::filesystem::remove(path + ".5");
    parallel_chunk_items = usual_chunk;
}

/// Membership of all query k-mers must agree with the definition, for any mask.
TEST(FMS_INDEX2, MEMBERSHIP) {
    std::mt19937_64 rng(4);
    for (int round = 0; round < 60; ++round) {
        size_t k = std::vector<size_t>{2, 3, 5, 11, 13, 31}[rng() % 6];
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 4000), k, round % 2);
        auto represented = represented_kmers(ms, k);
        queries q = make_queries(rng, ms, k);

        for (auto mask_type : {fmsi2::mask_kind::plain, fmsi2::mask_kind::packed, fmsi2::mask_kind::rrr}) {
            // Streaming with and without the kLCP array.
            mbwt m = compute_mbwt<uint64_t>(ms, k, (round / 2 + (int)mask_type) % 2);
            fmsi2::index index;
            fmsi2::build(index, m, k, mask_type, 64 << (round % 3));
            for (unsigned table : {0, 2, 10}) {
                fmsi2::build_prefix_table(index, table);
                for (unsigned batch : {1, 3, 32}) {
                    fmsi2::lookup_stats stats;
                    std::vector<int8_t> single(q.sequences.size(), 0), streamed(q.sequences.size(), 0);
                    fmsi2::lookup_batch<false, int8_t>(index, q.sequences.data(), q.offsets.data(), q.offsets.size(),
                                                       single.data(), batch, stats);
                    fmsi2::lookup_streaming<false, int8_t>(index, q.sequences.data(), q.segments.data(),
                                                           q.segments.size(), streamed.data(), batch, stats);
                    for (uint32_t offset : q.offsets) {
                        bool want = represented.count(canonical(q.sequences.substr(offset, k)));
                        ASSERT_EQ((bool)single[offset], want) << "k " << k << " round " << round;
                        ASSERT_EQ((bool)streamed[offset], want) << "k " << k << " round " << round << " (streaming)";
                    }
                }
            }
        }
    }
}

/// With one set occurrence per k-mer, the hashes are a bijection from the k-mers to [0, number of k-mers).
TEST(FMS_INDEX2, HASHES) {
    std::mt19937_64 rng(5);
    for (int round = 0; round < 40; ++round) {
        size_t k = std::vector<size_t>{3, 5, 13, 31}[rng() % 4];
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 4000), k, true);
        auto represented = represented_kmers(ms, k);
        queries q = make_queries(rng, ms, k);

        mbwt m = compute_mbwt<uint64_t>(ms, k, round / 9 % 2);
        fmsi2::index index;
        fmsi2::build(index, m, k, (fmsi2::mask_kind)(round % 3), 64 << (round / 3 % 3));
        fmsi2::build_prefix_table(index, round % 3 * 4);
        fmsi2::lookup_stats stats;
        std::vector<int64_t> single(q.sequences.size(), -1), streamed(q.sequences.size(), -1);
        fmsi2::lookup_batch<true, int64_t>(index, q.sequences.data(), q.offsets.data(), q.offsets.size(),
                                           single.data(), 7, stats);
        fmsi2::lookup_streaming<true, int64_t>(index, q.sequences.data(), q.segments.data(), q.segments.size(),
                                               streamed.data(), 7, stats);
        std::map<std::string, int64_t> hash_of;
        std::set<int64_t> hashes;
        for (uint32_t offset : q.offsets) {
            std::string kmer = canonical(q.sequences.substr(offset, k));
            ASSERT_EQ(single[offset], streamed[offset]) << "k " << k << " round " << round;
            if (!represented.count(kmer)) {
                ASSERT_EQ(single[offset], -1);
                continue;
            }
            ASSERT_GE(single[offset], 0);
            ASSERT_LT(single[offset], (int64_t)represented.size());
            auto inserted = hash_of.insert({kmer, single[offset]});
            ASSERT_EQ(inserted.first->second, single[offset]) << "a k-mer has two hashes";
            hashes.insert(single[offset]);
        }
        // The superstring itself is among the queries, so all k-mers were looked up.
        EXPECT_EQ(hash_of.size(), represented.size());
        EXPECT_EQ(hashes.size(), represented.size());
    }
}

// Access inverts the hashes: the k-mer of every hash is represented, has this hash, and the batched
// access returns the same k-mers for every batch size.
TEST(FMS_INDEX2, ACCESS) {
    std::mt19937_64 rng(6);
    for (int round = 0; round < 45; ++round) {
        size_t k = std::vector<size_t>{3, 5, 13, 31}[rng() % 4];
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 4000), k, true);
        auto represented = represented_kmers(ms, k);

        mbwt m = compute_mbwt<uint64_t>(ms, k, false);
        fmsi2::index index;
        fmsi2::build(index, m, k, (fmsi2::mask_kind)(round % 3), 64 << (round / 3 % 3));
        fmsi2::build_prefix_table(index, round % 3 * 4);
        ASSERT_EQ(index.mask_ones, represented.size()) << "k " << k << " round " << round;
        if (round % 2) fmsi2::build_access_support(index);

        std::vector<uint64_t> hashes(index.mask_ones);
        for (size_t h = 0; h < hashes.size(); ++h) hashes[h] = h;
        std::shuffle(hashes.begin(), hashes.end(), rng);
        std::string single(hashes.size() * k, 0);
        for (size_t i = 0; i < hashes.size(); ++i) fmsi2::access(index, hashes[i], single.data() + i * k);
        for (unsigned batch : {1u, 2u, 7u, 32u, 1000u}) {
            std::string batched(hashes.size() * k, 0);
            fmsi2::access_batch(index, hashes.data(), hashes.size(), batched.data(), batch);
            ASSERT_EQ(batched, single) << "k " << k << " round " << round << " batch " << batch;
        }

        // Looking the k-mers up gives the hashes back.
        std::vector<uint32_t> offsets(hashes.size());
        for (size_t i = 0; i < hashes.size(); ++i) offsets[i] = i * k;
        std::vector<int64_t> found(single.size() + 1, -1);
        fmsi2::lookup_stats stats;
        fmsi2::lookup_batch<true, int64_t>(index, single.data(), offsets.data(), offsets.size(), found.data(), 7, stats);
        std::set<std::string> seen;
        for (size_t i = 0; i < hashes.size(); ++i) {
            std::string kmer = single.substr(i * k, k);
            ASSERT_TRUE(represented.count(canonical(kmer))) << "k " << k << " round " << round << " k-mer " << kmer;
            ASSERT_EQ(found[offsets[i]], (int64_t)hashes[i]) << "k " << k << " round " << round << " k-mer " << kmer;
            seen.insert(canonical(kmer));
        }
        EXPECT_EQ(seen.size(), represented.size());
    }
}

/// The MBWT of the doubled text (doubled.h) computed naively: the suffixes of X 0 Y sorted with the end of
/// the text before the separator and the separator before the characters.
mbwt naive_doubled_mbwt(const std::string &ms, size_t k) {
    const size_t n = ms.size();
    std::string upper = ms;
    for (auto &c : upper) c = toupper(c);
    // '!' is the separator; the end of the text sorts before it as the end of a string.
    std::string text = upper + "!" + reverse_complement(upper);
    auto rank = [](char c) { return c == '!' ? 0 : 1 + (int)nucleotideToInt[(uint8_t)c]; };
    std::vector<size_t> sa(text.size() + 1);
    for (size_t i = 0; i < sa.size(); ++i) sa[i] = i;
    std::sort(sa.begin(), sa.end(), [&](size_t a, size_t b) {
        for (;; ++a, ++b) {
            if (b == text.size()) return false;
            if (a == text.size()) return true;
            if (text[a] != text[b]) return rank(text[a]) < rank(text[b]);
        }
    });
    auto represented = [&](size_t p) {
        if (p < n) return is_upper(ms[p]);
        if (p == n || p == text.size()) return false;
        size_t j = p - n - 1;
        return j + k <= n && is_upper(ms[n - j - k]);
    };
    auto has_prefix = [&](size_t p) { return p + k - 1 <= n || (p > n && p + k - 1 <= text.size()); };
    mbwt m;
    m.bwt.resize(sa.size());
    m.sa_transformed_mask = sdsl::bit_vector(sa.size(), 0);
    m.klcp = sdsl::bit_vector(sa.size(), 0);
    for (size_t i = 0; i < sa.size(); ++i) {
        size_t p = sa[i];
        if (p == 0) {
            m.dollar_position = i;
        } else if (p == n + 1) {
            m.dollar_position2 = i;
        } else {
            m.bwt[i] = nucleotideToInt[(uint8_t)text[p - 1]];
        }
        m.sa_transformed_mask[i] = represented(p);
        if (i + 1 < sa.size() && has_prefix(p) && has_prefix(sa[i + 1]) &&
            text.compare(p, k - 1, text, sa[i + 1], k - 1) == 0) {
            m.klcp[i] = 1;
        }
    }
    return m;
}

void expect_same_mbwt(const mbwt &a, const mbwt &b, const std::string &where) {
    ASSERT_EQ(a.bwt, b.bwt) << where;
    ASSERT_EQ(a.sa_transformed_mask, b.sa_transformed_mask) << where;
    ASSERT_EQ(a.klcp, b.klcp) << where;
    ASSERT_EQ(a.dollar_position, b.dollar_position) << where;
    ASSERT_EQ(a.dollar_position2, b.dollar_position2) << where;
}

/// The doubled text: its MBWT is the one of the suffixes sorted naively, and does not depend on the threads.
TEST(FMS_INDEX2, DOUBLED_MATCHES_NAIVE) {
    std::mt19937_64 rng(7);
    size_t usual_chunk = parallel_chunk_items;
    parallel_chunk_items = 64;
    for (int round = 0; round < 40; ++round) {
        size_t k = 2 + rng() % 31;
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % (round < 20 ? 300 : 3000)), k, round % 2);
        mbwt naive = naive_doubled_mbwt(ms, k);
        for (unsigned threads : {1u, 3u}) {
            mbwt m = doubled::compute_mbwt(ms, k, true, threads);
            expect_same_mbwt(m, naive, "k " + std::to_string(k) + " round " + std::to_string(round));
        }
    }
    parallel_chunk_items = usual_chunk;
}

/// The lean construction of the doubled text (merging the suffix arrays of X and Y) gives the same MBWT.
TEST(FMS_INDEX2, DOUBLED_LEAN_SAME_AS_LIBSAIS) {
    std::mt19937_64 rng(10);
    size_t usual_chunk = parallel_chunk_items;
    std::string temporary = (std::filesystem::temp_directory_path() / "fmsi2_lean_test").string();
    for (int round = 0; round < 60; ++round) {
        size_t k = 2 + rng() % 31;
        size_t length = round < 10 ? 1 + rng() % 5 : k + rng() % (round < 30 ? 300 : 20000);
        std::string ms = random_mask(rng, random_superstring(rng, length), k, round % 2);
        // Repeats of a few letters, where many suffixes of Y fall between the same suffixes of X.
        if (round % 7 == 3) ms = random_mask(rng, std::string(length / 2, 'A') + std::string(length - length / 2, 'C'), k, false);
        bool klcp = round % 3 != 0;
        mbwt want = doubled::compute_mbwt(ms, k, klcp, 1);
        for (unsigned threads : {1u, 3u, 8u}) {
            parallel_chunk_items = threads == 1 ? usual_chunk : 64;
            mbwt lean = doubled::compute_mbwt_lean(ms, k, klcp, threads, temporary);
            expect_same_mbwt(lean, want, "k " + std::to_string(k) + " round " + std::to_string(round) + " threads " +
                                             std::to_string(threads));
        }
        parallel_chunk_items = usual_chunk;
        ASSERT_FALSE(std::filesystem::exists(temporary + ".sa_x"));
    }
}

/// With both strands indexed, membership agrees with the definition and with the index of one strand.
TEST(FMS_INDEX2, FAST_MEMBERSHIP) {
    std::mt19937_64 rng(8);
    for (int round = 0; round < 40; ++round) {
        size_t k = std::vector<size_t>{2, 3, 5, 11, 13, 31}[rng() % 6];
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 4000), k, round % 2);
        auto represented = represented_kmers(ms, k);
        queries q = make_queries(rng, ms, k);
        for (auto mask_type : {fmsi2::mask_kind::plain, fmsi2::mask_kind::packed, fmsi2::mask_kind::rrr}) {
            mbwt m = doubled::compute_mbwt(ms, k, (round / 2 + (int)mask_type) % 2);
            fmsi2::index index;
            fmsi2::build(index, m, k, mask_type, 64 << (round % 3));
            ASSERT_EQ(index.strands, 2u);
            for (unsigned table : {0, 2, 10}) {
                fmsi2::build_prefix_table(index, table);
                for (unsigned batch : {1, 3, 32}) {
                    fmsi2::lookup_stats stats;
                    std::vector<int8_t> single(q.sequences.size(), 0), streamed(q.sequences.size(), 0);
                    fmsi2::lookup_batch<false, int8_t>(index, q.sequences.data(), q.offsets.data(), q.offsets.size(),
                                                       single.data(), batch, stats);
                    fmsi2::lookup_streaming<false, int8_t>(index, q.sequences.data(), q.segments.data(),
                                                           q.segments.size(), streamed.data(), batch, stats);
                    for (uint32_t offset : q.offsets) {
                        bool want = represented.count(canonical(q.sequences.substr(offset, k)));
                        ASSERT_EQ((bool)single[offset], want) << "k " << k << " round " << round;
                        ASSERT_EQ((bool)streamed[offset], want) << "k " << k << " round " << round << " (streaming)";
                    }
                }
            }
        }
    }
}

/// With both strands indexed and one set occurrence per k-mer, the hashes are a bijection from the represented
/// k-mers in both orientations to [0, 2 * number of k-mers), and access inverts them.
TEST(FMS_INDEX2, FAST_HASHES_AND_ACCESS) {
    std::mt19937_64 rng(9);
    for (int round = 0; round < 30; ++round) {
        size_t k = std::vector<size_t>{3, 5, 13, 31}[rng() % 4];
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 4000), k, true);
        auto represented = represented_kmers(ms, k);
        queries q = make_queries(rng, ms, k);

        mbwt m = doubled::compute_mbwt(ms, k, round / 9 % 2);
        fmsi2::index index;
        fmsi2::build(index, m, k, (fmsi2::mask_kind)(round % 3), 64 << (round / 3 % 3));
        fmsi2::build_prefix_table(index, round % 3 * 4);
        ASSERT_EQ(index.mask_ones, 2 * represented.size()) << "k " << k << " round " << round;
        fmsi2::lookup_stats stats;
        std::vector<int64_t> single(q.sequences.size(), -1), streamed(q.sequences.size(), -1);
        fmsi2::lookup_batch<true, int64_t>(index, q.sequences.data(), q.offsets.data(), q.offsets.size(),
                                           single.data(), 7, stats);
        fmsi2::lookup_streaming<true, int64_t>(index, q.sequences.data(), q.segments.data(), q.segments.size(),
                                               streamed.data(), 7, stats);
        std::map<std::string, int64_t> hash_of;
        std::set<int64_t> hashes;
        for (uint32_t offset : q.offsets) {
            std::string kmer = q.sequences.substr(offset, k);
            ASSERT_EQ(single[offset], streamed[offset]) << "k " << k << " round " << round;
            if (!represented.count(canonical(kmer))) {
                ASSERT_EQ(single[offset], -1);
                continue;
            }
            ASSERT_GE(single[offset], 0);
            ASSERT_LT(single[offset], 2 * (int64_t)represented.size());
            auto inserted = hash_of.insert({kmer, single[offset]});
            ASSERT_EQ(inserted.first->second, single[offset]) << "a k-mer has two hashes";
            hashes.insert(single[offset]);
            // Access gives the k-mer back, as it is.
            std::string accessed(k, 0);
            fmsi2::access(index, single[offset], accessed.data());
            ASSERT_EQ(accessed, kmer) << "k " << k << " round " << round;
        }
        // The superstring and its reverse complement are among the queries, so all k-mers were looked up.
        EXPECT_EQ(hash_of.size(), 2 * represented.size());
        EXPECT_EQ(hashes.size(), 2 * represented.size());
        if (round % 2) fmsi2::build_access_support(index);
        std::vector<uint64_t> all(index.mask_ones);
        for (size_t h = 0; h < all.size(); ++h) all[h] = h;
        std::string batched(all.size() * k, 0);
        fmsi2::access_batch(index, all.data(), all.size(), batched.data(), 5);
        for (size_t h = 0; h < all.size(); ++h) {
            std::string kmer = batched.substr(h * k, k);
            ASSERT_TRUE(hash_of.count(kmer)) << "k " << k << " round " << round;
            ASSERT_EQ(hash_of[kmer], (int64_t)h);
        }
    }
}

/// Streaming on an index of both strands (a loop of its own for k <= 32) agrees with the lookups on mosaic
/// reads: pieces of the superstring on both strands, random stretches and substitutions at several rates, more
/// reads than slots. The hashes agree too if every k-mer has one set occurrence.
TEST(FMS_INDEX2, FAST_STREAMING_MOSAICS) {
    std::mt19937_64 rng(11);
    for (int round = 0; round < 24; ++round) {
        size_t k = std::vector<size_t>{3, 5, 11, 15, 21, 31, 32}[rng() % 7];
        bool single = round % 2;
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 30000), k, single);
        std::string superstring = ms;
        for (auto &c : superstring) c = toupper(c);
        auto represented = represented_kmers(ms, k);
        queries q;
        for (int r = 0; r < 300; ++r) {
            size_t length = k + rng() % 300;
            unsigned per_mille = std::vector<unsigned>{0, 5, 20, 100}[rng() % 4];
            std::string read;
            while (read.size() < length) {
                size_t piece = 1 + rng() % 120;
                if (rng() % 4 == 0) {
                    read += random_dna(rng, piece);
                } else {
                    std::string s = superstring.substr(rng() % superstring.size(), piece);
                    read += rng() % 2 ? reverse_complement(s) : s;
                }
            }
            read.resize(length);
            for (auto &c : read) {
                if (rng() % 1000 < per_mille) c = "ACGT"[rng() % 4];
            }
            q.segments.push_back({(uint32_t)q.sequences.size(), (uint32_t)read.size()});
            for (size_t i = 0; i + k <= read.size(); ++i) q.offsets.push_back(q.sequences.size() + i);
            q.sequences += read;
        }
        for (bool klcp : {false, true}) {
            mbwt m = doubled::compute_mbwt(ms, k, klcp);
            fmsi2::index index;
            fmsi2::build(index, m, k, (fmsi2::mask_kind)(round / 2 % 3), 64 << (round % 3));
            for (unsigned table : {0u, 4u, 12u}) {
                fmsi2::build_prefix_table(index, table);
                for (unsigned batch : {1u, 5u, 64u}) {
                    fmsi2::lookup_stats stats;
                    std::vector<int64_t> looked_up(q.sequences.size(), -1), streamed(q.sequences.size(), -1);
                    fmsi2::lookup_batch<true, int64_t>(index, q.sequences.data(), q.offsets.data(), q.offsets.size(),
                                                       looked_up.data(), batch, stats);
                    fmsi2::lookup_streaming<true, int64_t>(index, q.sequences.data(), q.segments.data(),
                                                           q.segments.size(), streamed.data(), batch, stats);
                    ASSERT_EQ(stats.kmers, 2 * q.offsets.size());
                    for (uint32_t offset : q.offsets) {
                        bool want = represented.count(canonical(q.sequences.substr(offset, k)));
                        std::string where = "k " + std::to_string(k) + " round " + std::to_string(round) + " kLCP " +
                                            std::to_string(klcp) + " table " + std::to_string(table) + " batch " +
                                            std::to_string(batch) + " offset " + std::to_string(offset);
                        ASSERT_EQ(looked_up[offset] >= 0, want) << where;
                        ASSERT_EQ(streamed[offset] >= 0, want) << where << " (streaming)";
                        if (single) ASSERT_EQ(streamed[offset], looked_up[offset]) << where << " (hash)";
                    }
                }
            }
        }
    }
}

TEST(FMS_INDEX2, SAVE_AND_LOAD) {
    std::mt19937_64 rng(6);
    size_t k = 13;
    std::string ms = random_mask(rng, random_superstring(rng, 3000), k, true);
    queries q = make_queries(rng, ms, k);
    for (int variant = 0; variant < 6; ++variant) {
        auto mask_type = (fmsi2::mask_kind)(variant % 3);
        mbwt m = variant < 3 ? compute_mbwt<uint64_t>(ms, k, false) : doubled::compute_mbwt(ms, k, false);
        fmsi2::index built, loaded;
        fmsi2::build(built, m, k, mask_type, mask_type == fmsi2::mask_kind::plain ? 128 : 256);
        std::string path = testing::TempDir() + "fmsi2_test.fmsi2";
        fmsi2::save(built, path);
        fmsi2::load(loaded, path);
        std::filesystem::remove(path);
        EXPECT_EQ(built.strands, loaded.strands);
        EXPECT_EQ(built.dollar_position2, loaded.dollar_position2);
        fmsi2::build_prefix_table(built, 6);
        fmsi2::build_prefix_table(loaded, 6);
        fmsi2::lookup_stats stats;
        std::vector<int64_t> want(q.sequences.size(), -1), got(q.sequences.size(), -1);
        fmsi2::lookup_batch<true, int64_t>(built, q.sequences.data(), q.offsets.data(), q.offsets.size(), want.data(),
                                           16, stats);
        fmsi2::lookup_batch<true, int64_t>(loaded, q.sequences.data(), q.offsets.data(), q.offsets.size(), got.data(),
                                           16, stats);
        EXPECT_EQ(want, got);
    }
}

/// The lookups of an index: the hashes of all query k-mers.
std::vector<int64_t> hashes_of(fmsi2::index &index, const queries &q) {
    fmsi2::build_prefix_table(index, 6);
    fmsi2::lookup_stats stats;
    std::vector<int64_t> result(q.sequences.size(), -1);
    fmsi2::lookup_batch<true, int64_t>(index, q.sequences.data(), q.offsets.data(), q.offsets.size(), result.data(), 16, stats);
    return result;
}

/// An index file gives the same index in every layout, whichever way it was written, with its kLCP array.
TEST(FMS_INDEX2, FILE_LAYOUTS) {
    std::mt19937_64 rng(7);
    const std::string path = testing::TempDir() + "fmsi2_layouts_test.fmsi";
    for (int round = 0; round < 8; ++round) {
        size_t k = std::vector<size_t>{3, 13, 31, 33}[round % 4];
        std::string ms = random_mask(rng, random_superstring(rng, k + 100 + rng() % 70000), k, round % 2);
        queries q = make_queries(rng, ms, k);
        const bool both = round >= 4;
        mbwt_with_same_kmer = true;
        mbwt m = both ? doubled::compute_mbwt(ms, k, true) : k <= 32 ? compute_mbwt<uint64_t>(ms, k, true) : compute_mbwt<__uint128_t>(ms, k, true);
        mbwt_with_same_kmer = false;
        const sdsl::bit_vector klcp = m.klcp;
        const fmsi2::mask_statistics statistics = fmsi2::classify_mask(m);
        const fmsi2::description d = fmsi2::save_mbwt(m, k, both ? "fast" : "small", fmsi2::preset_layout(both ? "fast" : "small"),
                                                      statistics, path, 1 + round % 3);
        fmsi2::save_klcp_bits(klcp, k, d.uid(), path + ".klcp", 1 + round % 3);
        fmsi2::index built;
        fmsi2::build(built, m, k, fmsi2::mask_kind::plain);
        const std::vector<int64_t> want = hashes_of(built, q);

        for (uint32_t block_bytes : {0, 64, 128, 256}) {
            for (int mask : {-1, 0, 1, 2}) {
                fmsi2::layout layout;
                layout.bwt_block_bytes = block_bytes;
                layout.mask = mask;
                fmsi2::index loaded;
                fmsi2::load(loaded, path, 1 + (block_bytes / 64) % 3, layout);
                ASSERT_EQ(loaded.k, k);
                ASSERT_EQ(loaded.strands, both ? 2u : 1u);
                ASSERT_EQ(loaded.size, built.size);
                ASSERT_EQ(loaded.mask_ones, built.mask_ones);
                ASSERT_EQ(loaded.preset, both ? "fast" : "small");
                ASSERT_EQ(loaded.mask_type_name, statistics.type);
                // The recommended layout is that of the preset.
                if (block_bytes == 0) ASSERT_EQ(loaded.bwt.block_bytes(), both ? 64u : 256u);
                if (mask < 0) ASSERT_EQ(loaded.mask_type, both ? fmsi2::mask_kind::plain : fmsi2::mask_kind::packed);
                for (size_t i = 0; i < loaded.size; i += 1 + rng() % 7) ASSERT_EQ(loaded.mask_get(i), built.mask_get(i)) << "row " << i;
                ASSERT_EQ(hashes_of(loaded, q), want) << "round " << round << " blocks " << block_bytes << " mask " << mask;
                fmsi2::load_klcp(loaded, path + ".klcp", 2);
                ASSERT_EQ(loaded.klcp, klcp);

                // An index written from memory is the same file content-wise: it loads to the same index.
                if (block_bytes == 128 && mask == (int)fmsi2::mask_kind::packed) {
                    fmsi2::save(loaded, path + ".again", 2);
                    fmsi2::index again;
                    fmsi2::load(again, path + ".again");
                    ASSERT_EQ(again.bwt.block_bytes(), 128u);
                    ASSERT_EQ(again.mask_type, fmsi2::mask_kind::packed);
                    ASSERT_EQ(hashes_of(again, q), want);
                    // Its kLCP file is recognized by the identifier of the content.
                    fmsi2::load_klcp(again, path + ".klcp");
                    std::filesystem::remove(path + ".again");
                }
            }
        }
    }
    std::filesystem::remove(path);
    std::filesystem::remove(path + ".klcp");
}

/// Sections that a reader does not know are skipped unless they are required; damaged files are refused.
TEST(FMS_INDEX2, FILE_SECTIONS) {
    std::mt19937_64 rng(8);
    const size_t k = 11;
    std::string ms = random_mask(rng, random_superstring(rng, 5000000 / 100), k, true);
    queries q = make_queries(rng, ms, k);
    mbwt m = compute_mbwt<uint64_t>(ms, k, true);
    const sdsl::bit_vector klcp = m.klcp;
    const std::string path = testing::TempDir() + "fmsi2_sections_test.fmsi";
    const fmsi2::description d = fmsi2::save_mbwt(m, k, "small", fmsi2::preset_layout("small"), fmsi2::mask_statistics(), path);
    fmsi2::index reference;
    fmsi2::load(reference, path);
    const std::vector<int64_t> want = hashes_of(reference, q);
    ASSERT_EQ(reference.mask_type_name, "");

    // Copy the file section by section and add one that no reader knows.
    auto copy_with = [&](uint32_t flags) {
        fmsi2::file::reader in(path);
        fmsi2::file::writer out(path + ".more");
        for (const auto &e : in.sections()) {
            const std::string content = in.read_all(e);
            if (e.flags & fmsi2::file::DEFLATE) {
                out.add_stream(e.name(), e.encoding, e.flags, e.items, content.size(), true, 2,
                               [&](uint64_t begin, uint64_t end, uint8_t *buffer) { memcpy(buffer, content.data() + begin, end - begin); });
            } else {
                out.add(e.name(), e.encoding, e.flags, e.items, content.data(), content.size());
            }
        }
        const std::string future = "something of a later version";
        out.add("FUTURE", 7, flags, 1, future.data(), future.size());
        out.finish();
    };
    copy_with(0);
    {
        fmsi2::index loaded;
        fmsi2::load(loaded, path + ".more");
        ASSERT_EQ(hashes_of(loaded, q), want);
        fmsi2::file::reader in(path + ".more");
        ASSERT_EQ(in.sections().size(), 4u);
        ASSERT_EQ(in.read_all(in.get("FUTURE")), "something of a later version");
    }
    copy_with(fmsi2::file::REQUIRED);
    {
        fmsi2::index loaded;
        ASSERT_THROW(fmsi2::load(loaded, path + ".more"), std::invalid_argument);
    }
    std::filesystem::remove(path + ".more");

    // A kLCP file of another index is refused.
    {
        std::string other = random_mask(rng, random_superstring(rng, ms.size()), k, true);
        mbwt n = compute_mbwt<uint64_t>(other, k, true);
        const fmsi2::description e = fmsi2::save_mbwt(n, k, "small", fmsi2::preset_layout("small"), fmsi2::mask_statistics(), path + ".other");
        fmsi2::save_klcp_bits(n.klcp, k, e.uid(), path + ".other.klcp");
        fmsi2::save_klcp_bits(klcp, k, d.uid(), path + ".klcp");
        ASSERT_NE(d.uid(), e.uid());
        fmsi2::index loaded;
        fmsi2::load(loaded, path);
        ASSERT_THROW(fmsi2::load_klcp(loaded, path + ".other.klcp"), std::invalid_argument);
        fmsi2::load_klcp(loaded, path + ".klcp");
        ASSERT_EQ(loaded.klcp, klcp);
        std::filesystem::remove(path + ".other");
        std::filesystem::remove(path + ".other.klcp");
        std::filesystem::remove(path + ".klcp");
    }

    // Truncated and damaged files.
    std::ifstream in(path, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    auto write = [&](const std::string &content) {
        std::ofstream out(path + ".bad", std::ios::binary);
        out.write(content.data(), content.size());
    };
    for (size_t size : {(size_t)0, (size_t)10, (size_t)24, bytes.size() / 2, bytes.size() - 1}) {
        write(bytes.substr(0, size));
        fmsi2::index loaded;
        ASSERT_THROW(fmsi2::load(loaded, path + ".bad"), std::invalid_argument) << "size " << size;
    }
    {
        // A flipped bit of the compressed mask: the chunk does not decompress, or the mask does not match.
        fmsi2::file::reader table(path);
        const fmsi2::file::section_entry &mask = table.get("MASK");
        std::string damaged = bytes;
        damaged[mask.offset + mask.stored_bytes / 3] ^= 0x10;
        write(damaged);
        fmsi2::index loaded;
        ASSERT_THROW(fmsi2::load(loaded, path + ".bad"), std::invalid_argument);
        // A flipped bit of the BWT is found by the checksum.
        const fmsi2::file::section_entry &bwt = table.get("BWT");
        damaged = bytes;
        damaged[bwt.offset + bwt.stored_bytes / 2] ^= 0x01;
        write(damaged);
        fmsi2::file::reader bad(path + ".bad");
        ASSERT_THROW(bad.read(bad.get("BWT"), 1, [](uint64_t, const uint8_t *, size_t) {}, true), std::invalid_argument);
        bad.read(bad.get("MASK"), 1, [](uint64_t, const uint8_t *, size_t) {}, true);
    }
    std::filesystem::remove(path + ".bad");
    std::filesystem::remove(path);
}

/// The type of the mask that the constructions compute, against its definition.
TEST(FMS_INDEX2, MASK_TYPE) {
    std::mt19937_64 rng(9);
    auto computed = [&](const std::string &ms, size_t k, bool both) {
        mbwt_with_same_kmer = true;
        std::string copy = ms;
        mbwt m = both ? doubled::compute_mbwt(copy, k, false, 1 + rng() % 3) : compute_mbwt<uint64_t>(copy, k, false, false, 1 + rng() % 3);
        mbwt_with_same_kmer = false;
        return fmsi2::classify_mask(m);
    };
    // By the definition: the occurrences and the set occurrences of every k-mer, of both strands if [canonical].
    auto defined = [&](const std::string &ms, size_t k, bool canonical_kmers) {
        std::map<std::string, std::pair<size_t, size_t>> seen;
        for (size_t i = 0; i + k <= ms.size(); ++i) {
            std::string kmer = ms.substr(i, k);
            for (auto &c : kmer) c = toupper(c);
            auto &entry = seen[canonical_kmers ? canonical(kmer) : kmer];
            entry.first++;
            entry.second += is_upper(ms[i]);
        }
        bool single = true, uniform = true;
        size_t kmers = 0;
        for (const auto &item : seen) {
            if (item.second.second == 0) continue;
            ++kmers;
            single = single && item.second.second == 1;
            uniform = uniform && item.second.second == item.second.first;
        }
        return std::make_pair(std::string(single && uniform ? "min-one,max-one" : single ? "min-one" : uniform ? "max-one" : "other"), kmers);
    };
    std::map<std::string, int> seen_types;
    for (int round = 0; round < 200; ++round) {
        size_t k = std::vector<size_t>{2, 3, 5, 13, 31}[rng() % 5];
        std::string s = round % 5 == 4 ? random_dna(rng, k + rng() % 300) : random_superstring(rng, k + rng() % 600);
        std::string ms;
        switch (round % 4) {
        case 0:
            ms = random_mask(rng, s, k, true);
            break;
        case 1:
            ms = random_mask(rng, s, k, false);
            break;
        default: {
            // A max-one mask: a random set of canonical k-mers, all of whose occurrences are set.
            std::set<std::string> chosen;
            ms = s;
            for (auto &c : ms) c = tolower(c);
            for (size_t i = 0; i + k <= s.size(); ++i) {
                if (rng() % 2) chosen.insert(canonical(s.substr(i, k)));
            }
            for (size_t i = 0; i + k <= s.size(); ++i) {
                if (chosen.count(canonical(s.substr(i, k)))) ms[i] = toupper(ms[i]);
            }
            break;
        }
        }
        for (bool both : {false, true}) {
            const fmsi2::mask_statistics got = computed(ms, k, both);
            const auto want = defined(ms, k, both);
            ASSERT_EQ(got.source, both ? "computed" : "computed-one-strand");
            ASSERT_EQ(got.type, want.first) << ms << " k " << k << " both " << both;
            // A stated type is kept if the mask does not contradict it, and replaced with a warning if it does.
            for (const std::string stated : {"min-one", "max-one"}) {
                const fmsi2::mask_statistics by_header = fmsi2::decide_mask_type(got, "auto", stated);
                const fmsi2::mask_statistics by_user = fmsi2::decide_mask_type(got, stated, stated == "min-one" ? "max-one" : "");
                const bool holds = stated == "min-one" ? got.single : got.uniform;
                ASSERT_EQ(by_header.warning.empty(), holds);
                ASSERT_EQ(by_header.source, holds ? "header" : got.source);
                ASSERT_EQ(by_user.source, holds ? "declared" : got.source);
                const std::string kept = both && got.single && got.uniform ? "min-one,max-one" : stated;
                ASSERT_EQ(by_header.type, holds ? kept : got.type);
                ASSERT_EQ(by_user.type, holds ? kept : got.type);
            }
            if (!both) ASSERT_EQ(got.kmers, want.second) << ms;
            seen_types[got.type]++;
        }
    }
    // All four types occurred.
    ASSERT_EQ(seen_types.size(), 4u);
    // Without the rows of the same k-mers nothing is computed, and a stated type is taken as it is.
    mbwt none;
    ASSERT_EQ(fmsi2::classify_mask(none).type, "unknown");
    ASSERT_EQ(fmsi2::decide_mask_type(fmsi2::classify_mask(none), "max-one", "min-one").type, "max-one");
    ASSERT_EQ(fmsi2::decide_mask_type(fmsi2::classify_mask(none), "max-one", "min-one").source, "declared");
    ASSERT_EQ(fmsi2::decide_mask_type(fmsi2::classify_mask(none), "auto", "min-one").source, "header");
    ASSERT_EQ(fmsi2::decide_mask_type(fmsi2::classify_mask(none), "auto", "").type, "unknown");
    // Fields of a FASTA header, as KmerCamel writes them.
    const std::string header = "maskedsuperstring dataset='a b.fa' k=31 alg=greedy mask=min-one mode=bidirectional";
    ASSERT_EQ(header_field(header, "mask"), "min-one");
    ASSERT_EQ(header_field(header, "k"), "31");
    ASSERT_EQ(header_field(header, "mode"), "bidirectional");
    ASSERT_EQ(header_field(header, "ask"), "");
    ASSERT_EQ(header_field("k=13", "k"), "13");
}

/// The Tunstall-coded bits built from a stream of words, for densities from sparse to dense.
TEST(RANK2, PACKED_BITS_FROM_STREAM) {
    std::mt19937_64 rng(10);
    for (double density : {0.0, 0.001, 0.05, 0.3, 0.5, 0.77, 0.96, 0.999, 1.0}) {
        for (size_t size : {(size_t)0, (size_t)1, (size_t)63, (size_t)64, (size_t)65, (size_t)1000, (size_t)70001, (size_t)300000}) {
            std::vector<bool> bits(size);
            size_t ones = 0;
            for (size_t i = 0; i < size; ++i) {
                bits[i] = (rng() % 1000000) < density * 1000000;
                ones += bits[i];
            }
            fmsi2::tunstall_bits packed;
            // Groups of 1 to 64 bits.
            packed.build_stream(size, ones, [&](auto sink) {
                for (size_t i = 0; i < size;) {
                    const unsigned count = (unsigned)std::min<size_t>(size - i, i % 3 == 0 ? 64 : 1 + (i * 7) % 64);
                    uint64_t word = 0;
                    for (unsigned t = 0; t < count; ++t) word |= (uint64_t)bits[i + t] << t;
                    sink(word, count);
                    i += count;
                }
            });
            ASSERT_EQ(packed.size(), size);
            size_t rank = 0;
            for (size_t i = 0; i < size; ++i) {
                if (i % 17 == 0 || size < 2000) {
                    ASSERT_EQ(packed.rank(i), rank) << "density " << density << " size " << size << " position " << i;
                    ASSERT_EQ(packed[i], bits[i]) << "density " << density << " size " << size << " position " << i;
                }
                rank += bits[i];
            }
            ASSERT_EQ(packed.rank(size), ones);
        }
    }
}

/// The packed bytes of a sequence of 2-bit symbols, written and read in pieces.
TEST(RANK2, PACKED_SYMBOLS_IN_PIECES) {
    std::mt19937_64 rng(11);
    for (uint32_t block_bytes : {64, 128, 256}) {
        for (size_t size : {(size_t)1, (size_t)223, (size_t)224, (size_t)225, (size_t)100003}) {
            std::vector<uint8_t> symbols(size), bytes((size + 3) / 4, 0);
            for (size_t i = 0; i < size; ++i) {
                symbols[i] = rng() % 4;
                bytes[i / 4] |= symbols[i] << (2 * (i % 4));
            }
            fmsi2::quad_sequence sequence(size, block_bytes);
            for (size_t first = 0; first < bytes.size();) {
                const size_t count = std::min<size_t>(bytes.size() - first, 1 + rng() % 300);
                sequence.set_packed(first, bytes.data() + first, count);
                first += count;
            }
            sequence.finalize(2);
            size_t ranks[4] = {0, 0, 0, 0};
            for (size_t i = 0; i < size; ++i) {
                for (uint8_t c = 0; c < 4; ++c) ASSERT_EQ(sequence.rank(i, c), ranks[c]);
                ranks[symbols[i]]++;
            }
            std::vector<uint8_t> back(bytes.size(), 0xFF);
            for (size_t first = 0; first < bytes.size();) {
                const size_t count = std::min<size_t>(bytes.size() - first, 1 + rng() % 300);
                sequence.get_packed(first, back.data() + first, count);
                first += count;
            }
            ASSERT_EQ(back, bytes);
        }
    }
}

/// The code for texts beyond 2^31 characters and indexes beyond 2^32 rows (suffix arrays and tables of 64-bit
/// integers), reached with short texts by lowering the limits: the same MBWT, the same answers.
TEST(FMS_INDEX2, WIDE_INTEGERS) {
    std::mt19937_64 rng(12);
    for (int round = 0; round < 12; ++round) {
        size_t k = std::vector<size_t>{5, 13, 31, 33}[round % 4];
        std::string ms = random_mask(rng, random_superstring(rng, k + 200 + rng() % 30000), k, round % 2);
        queries q = make_queries(rng, ms, k);
        const bool both = round % 3 == 2;
        auto compute = [&]() {
            std::string copy = ms;
            return both ? doubled::compute_mbwt(copy, k, true, 1 + round % 3)
                        : k <= 32 ? compute_mbwt<uint64_t>(copy, k, true, false, 1 + round % 3)
                                  : compute_mbwt<__uint128_t>(copy, k, true, false, 1 + round % 3);
        };
        mbwt narrow = compute();
        suffix_array_32bit_limit = 0;
        mbwt wide = compute();
        suffix_array_32bit_limit = (size_t)INT32_MAX;
        ASSERT_EQ(narrow.bwt, wide.bwt);
        ASSERT_EQ(narrow.sa_transformed_mask, wide.sa_transformed_mask);
        ASSERT_EQ(narrow.klcp, wide.klcp);
        ASSERT_EQ(narrow.dollar_position, wide.dollar_position);
        ASSERT_EQ(narrow.dollar_position2, wide.dollar_position2);

        fmsi2::index a, b;
        fmsi2::build(a, narrow, k, fmsi2::mask_kind::plain);
        fmsi2::build(b, wide, k, fmsi2::mask_kind::plain);
        for (unsigned table : {0, 3, 7}) {
            fmsi2::build_prefix_table(a, table);
            fmsi2::prefix_table_32bit_limit = 0;
            fmsi2::build_prefix_table(b, table);
            fmsi2::prefix_table_32bit_limit = (uint64_t)1 << 32;
            ASSERT_TRUE(table == 0 || (!a.prefix_table32.empty() && a.prefix_table64.empty()));
            ASSERT_TRUE(table == 0 || (b.prefix_table32.empty() && !b.prefix_table64.empty()));
            fmsi2::lookup_stats stats;
            std::vector<int64_t> want(q.sequences.size(), -1), got(q.sequences.size(), -1), streamed(q.sequences.size(), -1);
            fmsi2::lookup_batch<true, int64_t>(a, q.sequences.data(), q.offsets.data(), q.offsets.size(), want.data(), 16, stats);
            fmsi2::lookup_batch<true, int64_t>(b, q.sequences.data(), q.offsets.data(), q.offsets.size(), got.data(), 16, stats);
            ASSERT_EQ(want, got) << "round " << round << " table " << table;
            std::vector<int8_t> member(q.sequences.size(), 0), member_streamed(q.sequences.size(), 0);
            fmsi2::lookup_batch<false, int8_t>(b, q.sequences.data(), q.offsets.data(), q.offsets.size(), member.data(), 16, stats);
            fmsi2::lookup_streaming<false, int8_t>(b, q.sequences.data(), q.segments.data(), q.segments.size(), member_streamed.data(), 8, stats);
            ASSERT_EQ(member, member_streamed) << "round " << round << " table " << table;
        }
    }
}

} // namespace
