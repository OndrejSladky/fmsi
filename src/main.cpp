#include <sdsl/suffix_arrays.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "compute_masks.h"
#include "functions.h"
#include "index.h"
#include "mask.h"
#include "parser.h"
#include "query.h"
#include "version.h"

static int usage() {
  fprintf(stderr, "HERE WILL BE USAGE\n");
  return 1;
}

static int usage_index() {
  fprintf(stderr, "HERE WILL BE USAGE\n");
  return 1;
}

static int usage_query() {
  fprintf(stderr, "HERE WILL BE USAGE\n");
  return 1;
}

static int version() {
  std::cout << VERSION << std::endl;
  return 0;
}

int ms_index(int argc, char *argv[]) {
  int c;
  bool usage = 0;
  int k;
  while ((c = getopt(argc, argv, "k:h")) >= 0) {
    switch (c) {
    case 'k':
      k = atoi(optarg);
      break;
    case 'h':
      usage = 1;
      break;
    default:
      usage_index();
      return 1;
    }
  }
  if (usage) {
    usage_index();
    return 0;
  }
  if (optind + 1 > argc) {
    usage_index();
    return 1;
  }
  std::string fn = std::string(argv[optind]);
  std::string superstring_path = fn + ".sstr";
  std::string mask_path = fn + ".mask";
  std::string index_path = fn + ".fm9";
  auto ms = read_masked_superstring(fn);
  // If k is not set, infer it assuming the standard format of the mask.
  if (!k)
    k = infer_k(ms.mask);
  write_superstring(superstring_path, ms.superstring);
  // Construct and dump the BW-transformed mask.
  bw_mask_t bw_transformed_mask = construct_bw_transformed_mask(ms);
  mask_dump(mask_path, bw_transformed_mask);
  // Construct and dump the FM-index.
  fm_index_t fm_index;
  sdsl::construct(fm_index, superstring_path, 1);
  sdsl::store_to_file(fm_index, index_path);
  return 0;
}

int ms_query(int argc, char *argv[]) {
  bool usage = false;
  int c;
  std::function<bool(size_t, size_t)> f = nullptr;
  while ((c = getopt(argc, argv, "f:h")) >= 0) {
    switch (c) {
    case 'f':
      f = mask_function(optarg);
      break;
    case 'h':
      usage = true;
      break;
    default:
      usage_index();
      return 1;
    }
  }
  if (usage) {
    usage_index();
    return 0;
  }
  if (optind + 2 > argc) {
    usage_query();
    return 1;
  }
  std::string fn = argv[optind];
  std::string mask_path = fn + ".mask";
  std::string index_path = fn + ".fm9";
  std::string kmer = argv[optind + 1];
  // Load FM-index.
  fm_index_t fm_index;
  sdsl::load_from_file(fm_index, index_path);
  // Load mask.
  bw_mask_t mask = mask_restore(mask_path);
  bool found;
  if (f == nullptr) {
    found = query(fm_index, mask, kmer);
  } else {
    bw_mask_rank_t rank(&mask);
    found = query_f(fm_index, mask, rank, f, kmer);
  }
  if (found)
    std::cout << "FOUND" << std::endl;
  else
    std::cout << "NOT FOUND" << std::endl;
  return 0;
}

int main(int argc, char *argv[]) {
  int ret = 0;
  if (argc < 2) {
    usage();
    return 0;
  }
  if (strcmp(argv[1], "index") == 0)
    ret = ms_index(argc - 1, argv + 1);
  else if (strcmp(argv[1], "query") == 0)
    ret = ms_query(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-v") == 0)
    return version();
  else
    return usage();

  return ret;
}
