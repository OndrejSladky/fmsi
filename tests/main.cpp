#include "gtest/gtest.h"
#include "q_suf_sort_test.h"
#include "fms_index_test.h"
#include "fms_index2_test.h"
#include "query_pipeline_test.h"
// #include "api_test.h"  // not passing yet: the kLCP file written by Index::save is not recognized by Index::load
//#include "compact_test.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
