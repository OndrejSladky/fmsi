#pragma once

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
                sdsl::rrr_vector<>({0, 0, 0, 1, 0, 1, 1, 1}),
                std::vector<size_t>({1, 3, 4, 7}),
                3
        };
        ret.ac_gt_rank.set_vector(&ret.ac_gt);
        ret.ac_rank.set_vector(&ret.ac);
        ret.gt_rank.set_vector(&ret.gt);
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

    TEST(FMS_INDEX, QUERY) {
        auto index = get_dummy_index();
        struct test_case {
            std::string query;
            bool want_result;
        };
        std::vector<test_case> tests = {
                {"A", true},
                {"AG", false},
                {"CA", true},
                {"GGTA", true},
                {"ATGG", false},
                {"GA", false},
                {"GGG", false},
                {"CC", true},
        };

        for (auto t: tests) {
            auto got_result = query(index, t.query, nullptr);

            EXPECT_EQ(got_result, t.want_result);
        }
    }

    TEST (FMS_INDEX, CONSTRUCT) {
        std::string masked_superstring = "CaGGTag";
        fms_index index = construct(masked_superstring);
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
}
