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
      {"A", "T"},
      {"CGT", "ACG"},
      {"CGCG", "CGCG"},
      {"ACTG", "CAGT"},
      {"AAAAC", "GTTTT"},
  };

  for (auto t : tests) {
    auto got_result = reverse_complement(t.input);

    EXPECT_EQ(got_result, t.want_result);
  }
}
} // namespace