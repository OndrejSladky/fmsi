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
            for (size_t distance : {0, 1, 2, 30, 250, 1100}) {
                if (i + distance > size) continue;
                size_t rank_i, rank_j;
                uint8_t c = rng() % 4;
                ranks.rank_pair(i, i + distance, c, rank_i, rank_j);
                ASSERT_EQ(rank_i, ranks.rank(i, c));
                ASSERT_EQ(rank_j, ranks.rank(i + distance, c));
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

} // namespace
