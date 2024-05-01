#include "gtest/gtest.h"
#include "fms_index_test.h"
#include "q_suf_sort_test.h"
//#include "compact_test.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
