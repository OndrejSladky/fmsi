#pragma once

#include <filesystem>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "../src/fmsi_api.h"
#include "gtest/gtest.h"

namespace {

/// The API against the definition of the represented k-mers: one k-mer, many k-mers, the k-mers of a sequence
/// (streamed or not), the dictionary, for both presets and all layouts, built in memory and loaded from a file.
TEST(FMSI_API, MEMBERSHIP_AND_DICTIONARY) {
    std::mt19937_64 rng(21);
    const std::string path = testing::TempDir() + "fmsi_api_test";
    for (int round = 0; round < 24; ++round) {
        const size_t k = std::vector<size_t>{3, 11, 31, 32, 33, 63}[round % 6];
        const std::string ms = random_mask(rng, random_superstring(rng, k + 50 + rng() % 5000), k, true);
        const auto represented = represented_kmers(ms, k);

        fmsi::BuildOptions options;
        options.preset = round % 2 ? fmsi::Preset::fast : round % 4 ? fmsi::Preset::extra_small : fmsi::Preset::small;
        options.layout = std::vector<fmsi::Layout>{fmsi::Layout::recommended, fmsi::Layout::compact, fmsi::Layout::fast}[round % 3];
        options.klcp = round % 4 < 2;
        options.access = round % 8 < 4;
        options.threads = 1 + round % 3;
        fmsi::Index built = fmsi::Index::build(ms, (int)k, options);
        const std::string file = built.save(path, 2);
        ASSERT_EQ(file, path + ".fmsi");
        fmsi::LoadOptions load;
        load.layout = std::vector<fmsi::Layout>{fmsi::Layout::fast, fmsi::Layout::recommended, fmsi::Layout::compact}[round % 3];
        load.access = round % 8 >= 4;
        fmsi::Index loaded = fmsi::Index::load(round % 2 ? file : path, load);

        for (const fmsi::Index *index : {&built, &loaded}) {
            ASSERT_EQ(index->k(), k);
            ASSERT_EQ(index->both_strands(), round % 2 == 1);
            ASSERT_EQ(index->has_klcp(), options.klcp);
            // A min-one mask: every k-mer of the set has one hash (two in an index of both strands, unless it is
            // its own reverse complement).
            ASSERT_TRUE(index->minimal_hashes());
            size_t palindromes = 0;
            for (const auto &kmer : represented) palindromes += kmer == reverse_complement(kmer);
            // (A k-mer that is its own reverse complement occurs twice among the rows of both strands.)
            ASSERT_EQ(index->number_of_hashes(), index->both_strands() ? 2 * represented.size() : represented.size());

            // A sequence with k-mers of the set, others, and characters that are not nucleotides.
            std::string sequence;
            for (int piece = 0; piece < 30; ++piece) {
                const size_t from = rng() % ms.size();
                std::string part = rng() % 3 ? ms.substr(from, 1 + rng() % (3 * k)) : random_dna(rng, 1 + rng() % (2 * k));
                for (auto &c : part) c = toupper(c);
                if (rng() % 2) part = reverse_complement(part);
                if (rng() % 5 == 0) part += "N";
                sequence += part;
            }
            const size_t count = sequence.size() >= k ? sequence.size() - k + 1 : 0;
            const std::vector<uint8_t> streamed = index->contains_all(sequence, true), independent = index->contains_all(sequence, false);
            const std::vector<int64_t> hashes = index->lookup_all(sequence, true), hashes_independent = index->lookup_all(sequence, false);
            ASSERT_EQ(streamed.size(), count);
            ASSERT_EQ(independent.size(), count);
            ASSERT_EQ(hashes.size(), count);
            std::string kmers;
            std::vector<bool> wanted;
            for (size_t i = 0; i < count; ++i) {
                const std::string kmer = sequence.substr(i, k);
                const bool want = kmer.find('N') == std::string::npos && represented.count(canonical(kmer));
                ASSERT_EQ((bool)streamed[i], want) << "round " << round << " k-mer " << kmer;
                ASSERT_EQ((bool)independent[i], want) << "round " << round << " k-mer " << kmer;
                ASSERT_EQ(hashes[i] >= 0, want);
                ASSERT_EQ(hashes[i], hashes_independent[i]);
                if (want) ASSERT_LT((uint64_t)hashes[i], index->number_of_hashes());
                if (i % 7 == 0) {
                    ASSERT_EQ(index->contains(kmer), want);
                    ASSERT_EQ(index->lookup(kmer), hashes[i]);
                    // The same k-mer in the other orientation: the same hash on one strand.
                    if (want && !index->both_strands()) ASSERT_EQ(index->lookup(reverse_complement(kmer)), hashes[i]);
                    if (want) ASSERT_EQ(canonical(index->access(hashes[i])), canonical(kmer));
                }
                kmers += kmer;
                wanted.push_back(want);
            }
            std::vector<uint8_t> present(count, 7);
            std::vector<int64_t> looked_up(count, 7);
            index->contains_batch(kmers.data(), count, present.data());
            index->lookup_batch(kmers.data(), count, looked_up.data());
            for (size_t i = 0; i < count; ++i) {
                ASSERT_EQ((bool)present[i], (bool)wanted[i]);
                ASSERT_EQ(looked_up[i], hashes[i]);
            }

            // The hashes are a bijection onto the k-mers of the set, and access inverts lookup.
            std::vector<uint64_t> all(index->number_of_hashes());
            for (size_t h = 0; h < all.size(); ++h) all[h] = h;
            std::string accessed(all.size() * k, '?');
            index->access_batch(all.data(), all.size(), accessed.data());
            std::set<std::string> seen;
            for (size_t h = 0; h < all.size(); ++h) {
                const std::string kmer = accessed.substr(h * k, k);
                ASSERT_TRUE(represented.count(canonical(kmer)));
                if (h % 5 == 0) ASSERT_EQ(index->access(h), kmer);
                if (index->both_strands()) {
                    // A k-mer that is its own reverse complement has two hashes, and lookup gives the smaller one.
                    const bool second = kmer == reverse_complement(kmer) && index->lookup(kmer) + 1 == (int64_t)h;
                    if (!second) ASSERT_EQ(index->lookup(kmer), (int64_t)h);
                } else {
                    ASSERT_EQ(index->lookup(kmer), (int64_t)h);
                    ASSERT_TRUE(seen.insert(canonical(kmer)).second);
                }
            }
            if (!index->both_strands()) ASSERT_EQ(seen.size(), represented.size());
            (void)palindromes;
        }
        ASSERT_THROW(built.contains(std::string(k + 1, 'A')), std::invalid_argument);
        ASSERT_THROW(built.access(built.number_of_hashes()), std::out_of_range);
    }
    std::filesystem::remove(path + ".fmsi");
    std::filesystem::remove(path + ".fmsi.klcp");
    ASSERT_THROW(fmsi::Index::load(path), std::invalid_argument);
}

} // namespace
