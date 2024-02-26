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
TEST(NORMALIZE, COUNT_K_MERS) {
  struct test_case {
    std::string superstring;
    mask_t mask;
    int k;
    assignable_function_t f;
    std::vector<std::string> want_result;
  };
  std::vector<test_case> tests = {
      {"ACGTAGATA", {1, 1,0, 0, 0, 1, 1, 0, 0}, 3, [](int ones, int total) {
         return ones % 2;
       },
       // Canonical $k$-mers in kmercount are the minimal from the back.
       {"ATA", "ATC"}
       }};

  for (auto &t : tests) {
      auto k_mers = camel::kh_init(S64);

      count_k_mers(k_mers, t.superstring, t.mask, t.k, t.f);
      auto k_mer_vec = camel::kMersToVec(k_mers);
      std::vector<std::string> k_mer_str_vec;
      for (auto k_mer: k_mer_vec) {
          k_mer_str_vec.push_back(camel::NumberToKMer(k_mer, t.k));
      }
      std::sort(k_mer_str_vec.begin(), k_mer_str_vec.end());

      EXPECT_EQ(k_mer_str_vec, t.want_result);
  }
}
} // namespace