#pragma once

#include <string>

#include "../src/parser.h"

#include "gtest/gtest.h"

namespace {
TEST(PARSER, REVERSE_COMPLEMENT) {
  struct test_case {
    std::string input;
    std::string want_result;
  };
  std::vector<test_case> tests = {
      {"A", "T"},       {"CGT", "ACG"},     {"CGCG", "CGCG"},
      {"ACTG", "CAGT"}, {"AAAAC", "GTTTT"},
  };

  for (auto t : tests) {
    auto got_result = reverse_complement(t.input);

    EXPECT_EQ(got_result, t.want_result);
  }
}
TEST(PARSER, COMPUTE_MASK_PATH) {
  struct test_case {
    std::string fn;
    int k;
    bool include_k;
    std::string want_result;
  };
  std::vector<test_case> tests = {
      {"human.fa", 31,  true,"human.fa.k31.mask"},
      {"bombyx.fa.xz", 8,  true, "bombyx.fa.xz.k8.mask"},
      {"bombyx.fa.xz", 8,  false, "bombyx.fa.xz.mask"},
  };

  for (auto t : tests) {
    auto got_result = compute_mask_path(t.fn, t.k, t.include_k);

    EXPECT_EQ(got_result, t.want_result);
  }
}
} // namespace