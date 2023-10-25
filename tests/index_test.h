#pragma once

#include <string>

#include "../src/index.h"

#include "gtest/gtest.h"

namespace {
TEST(INDEX, INDEX) {
  struct test_case {
    masked_superstring_t ms;
    mask_t want_result;
  };
  std::vector<test_case> tests = {
      // The input corresponds to mask 11010101100 (MS=CAcCtAgGTg$).
      {{{1, 1, 0, 1, 0, 1, 0, 1, 1, 0}, "CACCTAGGTG"},
       {0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1}}};

  for (auto t : tests) {
    auto got_result = construct_bw_transformed_mask(t.ms);

    for (size_t i = 0; i < t.want_result.size(); ++i) {
      EXPECT_EQ(t.want_result[i], got_result[i]);
    }
  }
}

TEST(INDEX, CONVERT_SUPERSTRING) {
  struct test_case {
    masked_superstring_t input;
    std::vector<qsint_t> want_result;
  };
  std::vector<test_case> tests = {{
      {mask_t{1, 0, 0, 1, 0, 1, 1, 0, 1, 1}, "CACCTAGGTG"},
      std::vector<qsint_t>({1, 0, 1, 1, 3, 0, 2, 2, 3, 2}),
  }};

  for (auto t : tests) {
    auto got_result = _convert_superstring(t.input);

    for (size_t i = 0; i < t.want_result.size(); ++i) {
      EXPECT_EQ(t.want_result[i], got_result[i]);
    }
  }
}
} // namespace