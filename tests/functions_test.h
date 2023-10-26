#pragma once

#include <string>

#include "../src/functions.h"

#include "gtest/gtest.h"

namespace {
TEST(FUNCTIONS, OR) {
  struct test_case {
    size_t ones;
    size_t total;
    bool want_result;
  };
  std::vector<test_case> tests = {
      {1, 5, true},
      {4, 4, true},
      {0, 3, false},
  };

  for (auto t : tests) {
    auto got_result = f_or(t.ones, t.total);

    EXPECT_EQ(got_result, t.want_result);
  }
}
TEST(FUNCTIONS, AND) {
  struct test_case {
    size_t ones;
    size_t total;
    bool want_result;
  };
  std::vector<test_case> tests = {
      {1, 5, false},
      {4, 4, true},
      {0, 3, false},
  };

  for (auto t : tests) {
    auto got_result = f_and(t.ones, t.total);

    EXPECT_EQ(got_result, t.want_result);
  }
}
TEST(FUNCTIONS, XOR) {
  struct test_case {
    size_t ones;
    size_t total;
    bool want_result;
  };
  std::vector<test_case> tests = {
      {1, 5, true},
      {4, 4, false},
      {0, 3, false},
  };

  for (auto t : tests) {
    auto got_result = f_xor(t.ones, t.total);

    EXPECT_EQ(got_result, t.want_result);
  }
}
TEST(FUNCTIONS, R_TO_S) {
  struct test_case {
    size_t ones;
    size_t total;
    size_t r;
    size_t s;
    bool want_result;
  };
  std::vector<test_case> tests = {
      {1, 5, 1, 3, true},  {4, 4, 3, 6, true},  {0, 3, 0, 0, true},
      {6, 7, 8, 9, false}, {4, 4, 2, 3, false}, {0, 3, 1, 15, false},
  };

  for (auto t : tests) {
    auto got_result = f_r_to_s(t.ones, t.total, t.r, t.s);

    EXPECT_EQ(got_result, t.want_result);
  }
}
TEST(FUNCTIONS, MASK_FUNCTION) {
  struct test_case {
    std::string name;
    size_t ones;
    size_t total;
    bool want_result;
  };
  std::vector<test_case> tests = {
      {"or", 1, 5, true},      {"and", 1, 5, false}, {"xor", 1, 5, true},
      {"1-7", 1, 5, true},     {"0-0", 1, 5, false}, {"1-1", 1, 5, true},
      {"30-110", 1, 5, false},
  };

  for (auto t : tests) {
    auto f = mask_function(t.name);
    auto got_result = f(t.ones, t.total);

    EXPECT_EQ(got_result, t.want_result);
  }
}
} // namespace