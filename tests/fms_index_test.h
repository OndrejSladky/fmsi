#pragma once

#include <sstream>

#include "../src/fms_index.h"

#include "gtest/gtest.h"

namespace {
    fms_index get_dummy_index() {
        fms_index ret = { //CAGGTAG$, 1011100$
                sdsl::bit_vector({1, 1, 0, 0, 0, 0, 1, 1}),
                sdsl::rank_support_v5<1>(),
                sdsl::bit_vector({1, 0, 0, 0}),
                sdsl::rank_support_v5<1>(),
                sdsl::bit_vector({0,1,0,0}),
                sdsl::rank_support_v5<1>(),
                sdsl::rrr_vector<255>({0, 0, 0, 1, 0, 1, 1, 1}),
                sdsl::rank_support_rrr<1, 255>(),
                std::vector<size_t>({1, 3, 4, 7}),
                3,
                sdsl::bit_vector({}),
        };
        ret.ac_gt_rank.set_vector(&ret.ac_gt);
        ret.ac_rank.set_vector(&ret.ac);
        ret.gt_rank.set_vector(&ret.gt);
        ret.mask_rank.set_vector(&ret.sa_transformed_mask);
        return ret;
    }
    fms_index get_dummy_index2() {
        fms_index ret = { //GGTAAGA$, 11001000$
                sdsl::bit_vector({0, 1, 1, 0, 0, 0, 1, 1, 1}),
                sdsl::rank_support_v5<1>(),
                sdsl::bit_vector({0, 0, 0, 0}),
                sdsl::rank_support_v5<1>(),
                sdsl::bit_vector({0,1,0,1,0}),
                sdsl::rank_support_v5<1>(),
                sdsl::rrr_vector<255>({0,0,1,0,0,1,1,0,0}),
                sdsl::rank_support_rrr<1, 255>(),
                std::vector<size_t>({1, 4, 4, 7}),
                5,
                sdsl::bit_vector({}),
        };
        ret.ac_gt_rank.set_vector(&ret.ac_gt);
        ret.ac_rank.set_vector(&ret.ac);
        ret.gt_rank.set_vector(&ret.gt);
        ret.mask_rank.set_vector(&ret.sa_transformed_mask);
        return ret;
    }
    fms_index get_dummy_index3() {
        fms_index ret = { // CACACAT$, 1110100$
                sdsl::bit_vector({1,0,0,0,0,0,0,0}),
                sdsl::rank_support_v5<1>(),
                sdsl::bit_vector({1,1,1,0,0,0,0}),
                sdsl::rank_support_v5<1>(),
                sdsl::bit_vector({1}),
                sdsl::rank_support_v5<1>(),
                sdsl::rrr_vector<255>({0,1,0,0,1,1,1,0}),
                sdsl::rank_support_rrr<1, 255>(),
                std::vector<size_t>({1, 4, 7, 7}),
                5,
                sdsl::bit_vector({0, 1, 0, 0, 1, 1, 0, 0}),
        };
        ret.ac_gt_rank.set_vector(&ret.ac_gt);
        ret.ac_rank.set_vector(&ret.ac);
        ret.gt_rank.set_vector(&ret.gt);
        ret.mask_rank.set_vector(&ret.sa_transformed_mask);
        return ret;
    }

    TEST(FMS_INDEX, RANK) {
        auto index = get_dummy_index();
        struct test_case {
            size_t i;
            byte c;
            size_t want_result;
        };
        std::vector<test_case> tests = {
                {4, 3, 1},
                {1, 3, 0},
                {1, 2, 1},
                {5, 0, 1},
                {7, 1, 1},
                {7, 2, 2},
                {8, 2, 3},
                {0, 2, 0},
        };

        for (auto t: tests) {
            auto got_result = rank(index, t.i, t.c);

            EXPECT_EQ(got_result, t.want_result);
        }
    }

    TEST(FMS_INDEX, RANK2) {
        auto index = get_dummy_index2();
        struct test_case {
            size_t i;
            byte c;
            size_t want_result;
        };
        std::vector<test_case> tests = {
                {4, 0, 2},
                {5, 0, 3},
                {6, 0, 3}
        };

        for (auto t: tests) {
            auto got_result = rank(index, t.i, t.c);

            EXPECT_EQ(got_result, t.want_result);
        }
    }

    TEST(FMS_INDEX, UPDATE_RANGE) {
        auto index = get_dummy_index();
        struct test_case {
            size_t i;
            size_t j;
            byte c;
            size_t want_i;
            size_t want_j;
        };
        std::vector<test_case> tests = {
                {0, 8, 0, 1, 3},
                {0, 5, 0, 1, 2},
                {4, 5, 0, 1, 2},
                {5, 6, 0, 2, 3},
                {0, 8, 1, 3, 4},
                {0, 8, 2, 4, 7},
                {0, 8, 3, 7, 8},
                {0, 2, 0, 1, 1},
        };

        for (auto t: tests) {
            update_range(index, t.i, t.j, t.c);

            EXPECT_EQ(t.i, t.want_i);
            EXPECT_EQ(t.j, t.want_j);
        }
    }

    TEST(FMS_INDEX, EXTEND_RANGE_WITH_KLCP) {
        auto index = get_dummy_index3();
        struct test_case {
            size_t i;
            size_t j;
            size_t want_i;
            size_t want_j;
        };
        std::vector<test_case> tests = {
                {4, 6, 4 ,7},
                {4, 5, 4, 7},
                {5, 6, 4, 7},
                {2, 3, 1, 3},
                {1, 2, 1, 3},
                {3, 4, 3, 4},
        };

        for (auto t: tests) {
            extend_range_with_klcp(index, t.i, t.j);

            EXPECT_EQ(t.i, t.want_i);
            EXPECT_EQ(t.j, t.want_j);
        }
    }

    TEST(FMS_INDEX, GET_RANGE_WITH_PATTERN) {
        auto index = get_dummy_index3();
        struct test_case {
            std::string pattern;
            size_t want_i;
            size_t want_j;
        };
        std::vector<test_case> tests = {
                {"ACA", 1, 3},
                {"CAC", 4, 6},
                {"CAT", 6, 7},
                {"AAA", 1, 1},
                {"TAC", 8, 8},
                {"A", 1, 4},
                {"CA", 4, 7},
                {"T", 7, 8},
        };

        for (auto t: tests) {
            auto pattern = (char*) t.pattern.data();
            size_t i, j;

            get_range_with_pattern(index, i, j, pattern, t.pattern.length());

            EXPECT_EQ(i, t.want_i);
            EXPECT_EQ(j, t.want_j);
        }
    }

    TEST(FMS_INDEX, KMER_ORDER_IF_PRESENT) {
        auto index = get_dummy_index3();
        struct test_case {
            size_t sa_start;
            size_t sa_end;
            int64_t want_result;
        };
        std::vector<test_case> tests = {
            {1, 2, 0},
            {2, 3, -1},
            {1, 1, -1},
            {3, 4, -1},
            {4, 6, 1},
            {5, 6, 2},
            {1, 6, 0},
        };
        for (auto t: tests) {

            int64_t got_result = kmer_order_if_present(index, t.sa_start, t.sa_end);

            EXPECT_EQ(got_result, t.want_result);
        }
    }


    TEST(FMS_INDEX, QUERY_KMERS_STREAMING) {
        auto index = get_dummy_index3();
        struct test_case {
            std::string query;
            int k;
            bool maximize_ones;
            std::string want_result;
        };
        std::vector<test_case> tests = {
                {"CACATACA",3, false, "111001"},
                {"TGTATGTG",3, false, "100111"},
                {"CACATTGT",3, false, "111001"},
                {"CACATACA",3, true,  "111001"}, // Assumes the first line is always selected to infer.
        };
        for (auto t: tests) {
            auto sequence = (char*) t.query.data();
            auto rc = ReverseComplementString(t.query.data(), t.query.length());
            std::stringstream got_result;

            if (t.maximize_ones)
                query_kmers_streaming<true>(index, sequence, rc, t.query.length(), t.k, false, got_result);
            else
                query_kmers_streaming<false>(index, sequence, rc, t.query.length(), t.k, false, got_result);

            EXPECT_EQ(got_result.str(), t.want_result);
        }
    }

    TEST(FMS_INDEX, QUERY_KMERS_STREAMING_ORDERS) {
        auto index = get_dummy_index3();
        struct test_case {
            std::string query;
            int k;
            bool maximize_ones;
            std::string want_result;
        };
        std::vector<test_case> tests = {
                {"CACATACA",3, false, "1,0,3,-1,-1,0"},
                {"TGTATGTG",3, false, "0,-1,-1,3,0,1"},
                {"CACATTGT",3, false, "1,0,3,-1,-1,0"},
        };
        for (auto t: tests) {
            auto sequence = (char*) t.query.data();
            auto rc = ReverseComplementString(t.query.data(), t.query.length());
            std::stringstream got_result;

            if (t.maximize_ones)
                query_kmers_streaming<true>(index, sequence, rc, t.query.length(), t.k, true, got_result);
            else
                query_kmers_streaming<false>(index, sequence, rc, t.query.length(), t.k, true, got_result);

            EXPECT_EQ(got_result.str(), t.want_result);
        }
    }

    TEST(FMS_INDEX, QUERY_ORDERS) {
        auto index = get_dummy_index();
        struct test_case {
            std::string query;
            int k;
            std::string want_result;
        };
        std::vector<test_case> tests = {
                {"A", 1, "3"},
                {"AG", 2, "-1"},
                {"CA", 2, "2"},
                {"GGTA", 4, "1"},
                {"ATGG", 4, "-1"},
                {"GA", 2, "-1"},
                {"GGG", 3, "-1"},
                {"CC", 2, "1"},
                {"CCAG", 2, "2,1,-1"},
        };

        for (auto t: tests) {
            std::stringstream got_result;
            
            query_kmers<query_mode::orr>(index, t.query.data(), t.query.length(), t.k, false, got_result, true);

            EXPECT_EQ(got_result.str(), t.want_result);
        }
    }

    TEST(FMS_INDEX, QUERY) {
        auto index = get_dummy_index();
        struct test_case {
            std::string query;
            std::string want_result;
        };
        std::vector<test_case> tests = {
                {"A", "1"},
                {"AG", "0"},
                {"CA", "1"},
                {"GGTA", "1"},
                {"ATGG", "0"},
                {"GA", "0"},
                {"GGG", "0"},
                {"CC", "1"},
        };

        for (auto t: tests) {
            std::stringstream got_result;
            
            query_kmers<query_mode::orr>(index, t.query.data(), t.query.length(), t.query.size(), false, got_result, false);

            EXPECT_EQ(got_result.str(), t.want_result);
        }
    }

    TEST(FMS_INDEX, QUERY2) {
        auto index = get_dummy_index2();
        struct test_case {
            std::string query;
            std::string want_result;
        };
        std::vector<test_case> tests = {
                {"AAGA", "1"},
                {"AAGAA", "0"},
                {"GGTTAAGA", "1"},
                {"GTTAAGA", "1"},
        };

        for (auto t: tests) {
            std::stringstream got_result;

            query_kmers<query_mode::orr>(index, t.query.data(), t.query.length(), t.query.size(), false, got_result, false);

            EXPECT_EQ(got_result.str(), t.want_result);
        }
    }

    TEST (FMS_INDEX, CONSTRUCT) {
        std::string masked_superstring = "CaGGTag";
        fms_index index = construct<int64_t>(masked_superstring, 31, false);
        fms_index want_index = get_dummy_index();
        EXPECT_EQ(index.sa_transformed_mask.size(), want_index.sa_transformed_mask.size());
        for (size_t i = 0; i < index.sa_transformed_mask.size(); ++i) {
            EXPECT_EQ(index.sa_transformed_mask[i], want_index.sa_transformed_mask[i]);
        }
        EXPECT_EQ(index.ac, want_index.ac);
        EXPECT_EQ(index.ac_gt, want_index.ac_gt);
        EXPECT_EQ(index.gt, want_index.gt);
        EXPECT_EQ(index.counts, want_index.counts);
        EXPECT_EQ(index.dollar_position, want_index.dollar_position);

    }

    TEST(FMS_INDEX, ACCESS) {
        auto index = get_dummy_index();
        struct test_case {
            size_t i;
            byte want_result;
        };
        std::vector<test_case> tests = {
                {0, 2},
                {1, 3},
                {2, 1},
                {3, 0},
                {4, 0},
                {5, 0},
                {6, 2},
                {7, 2},
        };

        for (auto t: tests) {
            auto got_result = access(index, t.i);

            EXPECT_EQ(got_result, t.want_result);
        }
    }

    TEST(FMS_INDEX, EXPORT_MS) {
        auto index = get_dummy_index();
        auto got_result = export_ms(index);
        std::string want_result = "CaGGTag";
        EXPECT_EQ(got_result, want_result);
    }

    TEST(FMS_INDEX, OBTAIN_KMER) {
        struct test_case {
            std::string masked_superstring;
            std::vector<int64_t> k_mers;
            size_t i;
            int64_t mask;
            int sparsity;
            int k;
            int64_t want_result;
        };
        std::vector<test_case> tests  = {
                {"CACACAT", {0b010001, 0b000100}, 0, 0b111111, 3, 3, 0b010001},
                {"CACACAT", {0b010001, 0b000100}, 1, 0b111111, 3, 3, 0b000100},
                {"CACACAT", {0b010001, 0b000100}, 2, 0b111111, 3, 3, 0b010001},
                {"CACACAT", {0b010001, 0b000100}, 3, 0b111111, 3, 3, 0b000100},
                {"CACACAT", {0b010001, 0b000100}, 4, 0b111111, 3, 3, 0b010011},
        };

        for (auto t: tests) {
            auto got_result = obtain_kmer(t.k_mers, t.masked_superstring, t.i, t.sparsity, t.mask, t.k);

            EXPECT_EQ(got_result, t.want_result);
        }
    }

    TEST(FMS_INDEX, CONSTRUCT_KLCP) {
        struct test_case {
            std::string masked_superstring;
            qsint_t* isa;
            size_t k;
            std::vector<bool> want_result;
        };
        std::vector<test_case> tests = {
                {
                        "CACACat", new qsint_t[8]{4, 1, 5, 2, 6, 3, 7, 0}, 3,
                        std::vector<bool> {0, 1, 0, 0, 1, 1, 0, 0},
                },
                {
                    "CACACat", new qsint_t[8]{4, 1, 5, 2, 6, 3, 7, 0}, 2 ,
                    std::vector<bool>{0, 1, 1, 0, 1, 1, 0, 0},
                },
                {
                    "AAAAT", new qsint_t[6]{1,2,3,4, 5, 0}, 3,
                    std::vector<bool>{0, 1, 1, 0, 0, 0},
                },
                {
                    "AAAAT", new qsint_t[6]{1,2,3,4, 5, 0}, 2,
                            std::vector<bool>{0, 1, 1, 1, 0, 0},
                },
        };

        for (auto t: tests) {
            qsint_t* sa = new qsint_t[t.masked_superstring.length() + 1];
            QSufSortGenerateSaFromInverse(t.isa, sa, (qsint_t)t.masked_superstring.length());

            auto got_result = construct_klcp<uint64_t>(sa, t.masked_superstring, t.k - 1);

            EXPECT_EQ(std::vector<bool>(got_result.begin(), got_result.end()), t.want_result);
        }

    }
}
