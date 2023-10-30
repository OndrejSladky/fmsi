#pragma once

#include "../src/compute_masks.h"

#include "gtest/gtest.h"

namespace {
TEST(COMPUTE_MASKS, COMPUTE_L_MASK) {
  struct test_case {
    mask_t input;
    int k;
    int l;
    mask_t want_result;
  };
  std::vector<test_case> tests = {
      {{1,0,1, 1, 0, 0}, 3, 3,
       {1,0,1, 1, 0, 0}},
      {{1,0,1, 1, 0, 0}, 3, 2,
              {1,1,1, 1, 1, 0}},
      {{1,0,1, 1, 0, 0}, 3, 1,
              {1,1,1, 1, 1, 1}},
      {{1,0,0, 1, 0, 0, 0}, 4, 3,
              {1,1,0, 1, 1, 0, 0}},
  };

  for (auto t : tests) {
    auto got_result = compute_l_mask(t.input, t.k, t.l);

    EXPECT_EQ(got_result, t.want_result);
  }
}
} // namespace