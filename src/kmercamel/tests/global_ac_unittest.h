#pragma once
#include "../src/ac/global_ac.h"

#include "gtest/gtest.h"

namespace {
    TEST(GlobalAC, Suffix) {
        EXPECT_EQ("CGT", Suffix(KMer{"ACGT"}, 1));
        EXPECT_EQ("CGTA", Suffix(KMer{"CGTA"}, 0));
        EXPECT_EQ("", Suffix(KMer{"CC"}, 2));
    }

    TEST(GlobalAC, SuperstringFromPath) {
        struct TestCase {
            std::vector<OverlapEdge> path;
            std::vector<KMer> kMers;
            int k;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        std::vector<OverlapEdge>{OverlapEdge{1, 0, 2},OverlapEdge{0, 2, 1}},
                        {KMer{"ACG"}, KMer{"TAC"}, KMer{"GGC"}},
                        3,
                        "TAcGgc",
                },
                {
                        std::vector<OverlapEdge>{OverlapEdge{2, 1, 2},OverlapEdge{1, 3, 1}},
                        {KMer{"GCC"}, KMer{"ACG"}, KMer{"TAC"}, KMer{"GGC"}, KMer{"CGT"}, KMer{"GTA"}},
                        3,
                        "TAcGgc",
                },
        };

        for (auto t : tests) {
            std::stringstream of;

            SuperstringFromPath(t.path, t.kMers, of, t.k);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(GlobalAC, OverlapHamiltonianPath) {
        struct TestCase {
            std::vector<KMer> kMers;
            std::vector<OverlapEdge> wantResult;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {
                        std::vector<KMer>{KMer{"ACG"}, KMer{"TAC"}, KMer{"GGC"}},
                        std::vector<OverlapEdge>{OverlapEdge{1, 0, 2},OverlapEdge{0, 2, 1}},
                        false,
                },
                {
                        {KMer{"TAA"}, KMer{"TTT"}, KMer{"TTA"}, KMer{"AAA"}},
                        {{1, 2, 2}, {0, 3, 2}},
                        true,
                },
        };

        for (auto t : tests) {
            std::vector<OverlapEdge> got = OverlapHamiltonianPathAC(t.kMers, t.complements);
            EXPECT_EQ(t.wantResult.size(), got.size());
            for (size_t i = 0; i < t.wantResult.size(); ++i) {
                EXPECT_EQ(t.wantResult[i].firstIndex, got[i].firstIndex);
                EXPECT_EQ(t.wantResult[i].secondIndex, got[i].secondIndex);
                EXPECT_EQ(t.wantResult[i].overlapLength, got[i].overlapLength);
            }
        }
    }


    TEST(GlobalAC, GlobalAC) {
        struct TestCase {
            std::string wantResult;
            std::vector<KMer> input;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {"TACgt", {KMer{"CGT"}, KMer{"TAC"}, KMer{"ACG"}}, false},
                {"GcTa", {KMer{"TA"}, KMer{"GC"}, }, false},
                {"ACgTtt", {KMer{"CGT"}, KMer{"TTT"}, KMer{"ACG"}}, false},
                {"TActt", {KMer{"TACT"}, KMer{"ACTT"}}, false},
                {"TActTaaGgac", {KMer{"TACT"}, KMer{"ACTT"}, KMer{"GGAC"}, KMer{"TAAG"}}, false},
                {"TTtcttttttttttttttttttttttttttga", {KMer{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}, KMer{"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"}}, false},
                {"AAcAaatCccc", {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}}, true},
        };

        for (auto &&t : tests) {
            std::stringstream of;

            GlobalAC(t.input, of,  t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}