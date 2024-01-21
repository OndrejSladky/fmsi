#pragma once

#include <string>

#include "../src/normalize.h"

#include "gtest/gtest.h"

namespace {
TEST(NORMALIZE, SEPARATE_MASK_AND_SUPERSTRING) {
  struct test_case {
    std::string input;
    masked_superstring_t want_result;
  };
  std::vector<test_case> tests = {
      {{"ACGtGTaa"}, {{1, 1, 1, 0, 1, 1, 0, 0}, "ACGTGTAA"}},
      {{"AgtaGaaa"}, {{1, 0, 0, 0, 1, 0, 0, 0}, "AGTAGAAA"}}};

  for (auto t : tests) {
    auto got_result = separate_mask_and_superstring(t.input);

    EXPECT_EQ(got_result.mask, t.want_result.mask);
  }
}
} // namespace