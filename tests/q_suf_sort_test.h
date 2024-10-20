#pragma once

#include "../src/QSufSort.h"

#include "gtest/gtest.h"

namespace {
TEST(EXTERN, QSUFSORT) {
  struct test_case {
    std::vector<qsint_t> input;
    std::vector<qsint_t> want_result;
  };
  std::vector<test_case> tests = {
      // The input is CACCTAGGTG($)
      // Its BWT is GCT$ACTAGCG
      {{1, 0, 1, 1, 3, 0, 2, 2, 3, 2},
       {3, 1, 4, 5, 9, 2, 7, 8, 10, 6, 0}}};

  for (auto t : tests) {

    auto workplace = (qsint_t *)malloc((t.input.size() + 1) * sizeof(qsint_t));
    auto input = (qsint_t *)malloc((t.input.size() + 1) * sizeof(qsint_t));
    for (size_t i = 0; i < t.input.size(); ++i) {
      input[i] = t.input[i];
    }

    QSufSortSuffixSort(input, workplace, (qsint_t)t.input.size(), 3, 0, 0);

    for (size_t i = 0; i < t.want_result.size(); ++i) {
      EXPECT_EQ(t.want_result[i], input[i]);
    }
  }
}

TEST(EXTERN, CONSTRUCT_SA) {
    struct test_case {
        std::vector<qsint_t> input;
        std::vector<qsint_t> want_result;
    };
    std::vector<test_case> tests = {
            {{3, 1, 4, 5, 9, 2, 7, 8, 10, 6, 0},
             {10, 1, 5, 0, 2, 3, 9, 6, 7, 4, 8}}
    };
    for (auto t : tests) {
        auto workplace = (qsint_t *)malloc((t.input.size() + 1) * sizeof(qsint_t));
        auto input = (qsint_t *)malloc((t.input.size() + 1) * sizeof(qsint_t));
        for (size_t i = 0; i < t.input.size(); ++i) {
            input[i] = t.input[i];
        }



        QSufSortGenerateSaFromInverse(input, workplace, (qsint_t)t.input.size() - 1);

        for (size_t i = 0; i < t.want_result.size(); ++i) {
            EXPECT_EQ(t.want_result[i], workplace[i]);
        }
    }
}

} // namespace
