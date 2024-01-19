#include "compute_masks_test.h"
#include "functions_test.h"
#include "index_test.h"
#include "parser_test.h"
#include "q_suf_sort_test.h"
#include "normalize_test.h"

#include "gtest/gtest.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}