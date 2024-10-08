#pragma once

#include "../src/lower_bound.h"

#include "gtest/gtest.h"

namespace {
    TEST(LowerBound, LowerBoundLength) {
        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            size_t wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        {KMerToNumber({"ACG"}), KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                        3, false, 5
                },
                {
                        {KMerToNumber({"AC"}), KMerToNumber({"GG"})},
                        2, true, 3
                },
                {
                        {KMerToNumber({"AAACCC"}), KMerToNumber({"CCCAAA"}), KMerToNumber({"TGGGGT"})},
                        6, false, 11
                },
        };

        for (auto &t : tests) {
            auto gotResult = LowerBoundLength(t.kMers, t.k, t.complements);

            ASSERT_EQ(t.wantResult, gotResult);
        }
    }
}