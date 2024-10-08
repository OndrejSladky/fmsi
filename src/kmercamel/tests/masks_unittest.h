#pragma once
#include "../src/masks.h"

#include <algorithm>
#include <vector>
#include <string>
#include <filesystem>

#include "gtest/gtest.h"

namespace {

// Retrieving current path on Windows does not work as on linux
// therefore the following unittest is linux-specific.
#ifdef __unix__

    TEST(Masks, OptimizeOnes) {
        std::string path = std::filesystem::current_path();
        path += "/tests/testdata/masktest.fa";

        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            bool minimize;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {
                    {KMerToNumber({"ACG"}), KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                    3, false, true,
                    "> superstring\nACgTaacgt\n"
                },
                {
                    {KMerToNumber({"ACG"}), KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                    3, false, false,
                    "> superstring\nACgTaACgt\n"
                },
                {
                    {KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                    3, true, true,
                    "> superstring\nAcgTaacgt\n"
                },
                {
                    {KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                    3, true, false,
                    "> superstring\nACgTaACgt\n"
                },
        };

        for (auto &t : tests) {
            std::stringstream of;
            std::ifstream in(path);
            if (!in.is_open()) {
                throw std::invalid_argument("couldn't open file " + path);
            }
            auto kMersDict = kh_init_S64();
            int ret;
            for (auto &kMer : t.kMers) kh_put_S64(kMersDict, kMer, &ret);

            OptimizeOnes(in, of, kMersDict, t.k, t.complements, t.minimize);
            in.close();

            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(Mask, OptimizeRuns) {
        std::string path = std::filesystem::current_path();

        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            bool approximate;
            std::string relativePath;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        {
                            KMerToNumber({"AC"}),
                            KMerToNumber({"CG"}),
                            KMerToNumber({"GT"}),
                            KMerToNumber({"TT"}),
                            KMerToNumber({"AT"})
                        },
                        2, false, false,
                        "/tests/testdata/runstest.fa",
                        "> superstring\nacgcgttACGtATt\n"
                },
                {
                        {
                                KMerToNumber({"AC"}),
                                KMerToNumber({"CG"}),
                                KMerToNumber({"GT"}),
                                KMerToNumber({"TT"}),
                                KMerToNumber({"AT"})
                        },
                        2, false, true,
                        "/tests/testdata/runstest.fa",
                        "> superstring\nACgCGTtACGtATt\n"
                },
                {
                        {
                                KMerToNumber({"ACT"}),
                                KMerToNumber({"CTA"}),
                                KMerToNumber({"GTA"})
                        },
                        3, true, false,
                        "/tests/testdata/runstest2.fa",
                        "> superstring\nACTAGta\n"
                },
        };

        for (auto &t : tests) {
            std::stringstream of;
            auto totalPath = path + t.relativePath;
            auto kMersDict = kh_init_S64();
            int ret;
            for (auto &kMer : t.kMers) kh_put_S64(kMersDict, kMer, &ret);

            OptimizeRuns(totalPath, kMersDict, of, t.k, t.complements, t.approximate);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
#endif
}
