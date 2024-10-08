#pragma once
#include "../src/global.h"

#include <algorithm>

#include "gtest/gtest.h"
typedef unsigned char byte;
namespace {
    TEST(Global, PartialPreSort) {
        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            std::vector<kmer_t> wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        {KMerToNumber({"GTA"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        3,
                        {KMerToNumber({"GGC"}), KMerToNumber({"GTA"}), KMerToNumber({"TAC"})},
                },
                {
                        {KMerToNumber({"TTTTTTTTTTTTT"}),KMerToNumber({"AAAAAAAAAAAAA"}),  KMerToNumber({"GCGCGCGCGCGCG"})},
                        13,
                        {KMerToNumber({"AAAAAAAAAAAAA"}), KMerToNumber({"GCGCGCGCGCGCG"}), KMerToNumber({"TTTTTTTTTTTTT"})},
                },
                {
                    // Sorting is done only based on the first 11 chars.
                        {KMerToNumber({"AAAAAAAAAAAAT"}),KMerToNumber({"AAAAAAAAAAAAA"})},
                        13,
                        {KMerToNumber({"AAAAAAAAAAAAT"}),KMerToNumber({"AAAAAAAAAAAAA"})},
                },
        };

        for (auto t : tests) {
            PartialPreSort(t.kMers, t.k);

            EXPECT_EQ(t.wantResult, t.kMers);
        }
    }

    TEST(Global, SuperstringFromPath) {
        struct TestCase {
            overlapPath path;
            std::vector<kmer_t> kMers;
            int k;
            std::string wantResult;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {
                        {{2, 0, (size_t)-1}, {1, 2, (byte)-1}},
                        std::vector<kmer_t>{KMerToNumber({"ACG"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        3,
                        "TAcGgc",
                        false,
                },
                {
                        {{4, 3, 1, (size_t)-1, 5, (size_t)-1}, {1 ,2, 1, (byte)-1, 2, (byte)-1}},
                        std::vector<kmer_t>{KMerToNumber({"GCC"}), KMerToNumber({"ACG"}), KMerToNumber({"TAC"})},
                        3,
                        "GcCGta",
                        true,
                },
        };

        for (auto t : tests) {
            std::stringstream of;

            SuperstringFromPath(t.path, t.kMers, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(Global, OverlapHamiltonianPath) {
        struct TestCase {
            std::vector<kmer_t> kMers;
            overlapPath wantResult;
            int k;
            bool complements;
            bool lower_bound;
        };
        std::vector<TestCase> tests = {
                {
                        {KMerToNumber({"AT"})},
                        {{(size_t)-1, (size_t)-1}, {(byte)-1, (byte)-1}},
                        2,
                        true,
                        false,
                },
                {
                        {KMerToNumber({"ACG"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        {{2, 0, (size_t)-1}, {1, 2, (byte)-1}},
                        3,
                        false,
                        false,
                },
                {
                        {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}), KMerToNumber({"AACA"})},
                        {{4, 3, 0, 5, (size_t)-1, (size_t)-1}, {2, 2, 3, 3, (byte)-1, (byte)-1}},
                        4,
                        true,
                        false,
                },
                {
                        {KMerToNumber({"ACG"}), KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                        {{1, 2, 0}, {2, 1, 1}},
                        3,
                        false,
                        true,
                },
                {
                        {KMerToNumber({"ACC"}), KMerToNumber({"CGG"})},
                        {{3, 2, 1, 0}, {2, 2, 0, 0}},
                        3,
                        true,
                        true,
                },
        };

        for (auto t : tests) {
            overlapPath got = OverlapHamiltonianPath( t.kMers, t.k, t.complements, t.lower_bound);
            EXPECT_EQ(t.wantResult.first, got.first);
            EXPECT_EQ(t.wantResult.second, got.second);
        }
    }

    TEST(Global, Global) {
        struct TestCase {
            std::string wantResult;
            int k;
            std::vector<kmer_t> input;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {"TACgt", 3, {KMerToNumber({"CGT"}), KMerToNumber({"TAC"}), KMerToNumber({"ACG"})}, false},
                {"ACgTtt", 3, {KMerToNumber({"CGT"}), KMerToNumber({"TTT"}), KMerToNumber({"ACG"})}, false},
                {"TActt", 4, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"})}, false},
                {"TActTaaGgac", 4, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"}), KMerToNumber({"GGAC"}), KMerToNumber({"TAAG"})}, false},
                {"TTtcttttttttttttttttttttttttttga", 31, {KMerToNumber({"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}), KMerToNumber({"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"})}, false},
                {"AtTTgtt", 4, {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}), KMerToNumber({"AACA"})}, true},
        };

        for (auto &&t : tests) {
            std::stringstream of;

            Global(t.input, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
