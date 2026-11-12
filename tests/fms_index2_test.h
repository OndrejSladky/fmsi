#pragma once

#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

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

TEST(RANK2, QUAD_RANK) {
    std::mt19937_64 rng(1);
    for (size_t size : {0, 1, 223, 224, 225, 4000, 57344, 57345, 130000}) {
        std::vector<uint8_t> symbols(size);
        fmsi2::quad_rank ranks(size);
        for (size_t i = 0; i < size; ++i) {
            symbols[i] = rng() % 4;
            ranks.set(i, symbols[i]);
        }
        ranks.finalize();
        size_t counts[4] = {0, 0, 0, 0};
        for (size_t i = 0; i <= size; ++i) {
            for (int c = 0; c < 4; ++c) ASSERT_EQ(ranks.rank(i, c), counts[c]) << "size " << size << " position " << i;
            for (size_t distance : {0, 1, 2, 30, 250}) {
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

/// Membership of all query k-mers must agree with the definition, for any mask.
TEST(FMS_INDEX2, MEMBERSHIP) {
    std::mt19937_64 rng(4);
    for (int round = 0; round < 60; ++round) {
        size_t k = std::vector<size_t>{2, 3, 5, 11, 13, 31}[rng() % 6];
        std::string ms = random_mask(rng, random_superstring(rng, k + rng() % 4000), k, round % 2);
        auto represented = represented_kmers(ms, k);
        queries q = make_queries(rng, ms, k);

        for (auto mask_type : {fmsi2::mask_kind::plain, fmsi2::mask_kind::rrr}) {
            mbwt m = compute_mbwt<uint64_t>(ms, k, false);
            fmsi2::index index;
            fmsi2::build(index, m, k, mask_type);
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

        mbwt m = compute_mbwt<uint64_t>(ms, k, false);
        fmsi2::index index;
        fmsi2::build(index, m, k, round % 2 ? fmsi2::mask_kind::plain : fmsi2::mask_kind::rrr);
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

TEST(FMS_INDEX2, SAVE_AND_LOAD) {
    std::mt19937_64 rng(6);
    size_t k = 13;
    std::string ms = random_mask(rng, random_superstring(rng, 3000), k, true);
    queries q = make_queries(rng, ms, k);
    for (auto mask_type : {fmsi2::mask_kind::plain, fmsi2::mask_kind::rrr}) {
        mbwt m = compute_mbwt<uint64_t>(ms, k, false);
        fmsi2::index built, loaded;
        fmsi2::build(built, m, k, mask_type);
        std::string path = testing::TempDir() + "fmsi2_test.fmsi2";
        fmsi2::save(built, path);
        fmsi2::load(loaded, path);
        std::filesystem::remove(path);
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
