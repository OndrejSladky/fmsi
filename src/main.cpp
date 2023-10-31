#include <sdsl/suffix_arrays.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <filesystem>

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

static int usage_clean() {
    fprintf(stderr, "HERE WILL BE USAGE\n");
    return 1;
}

static int version() {
  std::cout << VERSION << std::endl;
  return 0;
}

std::string compute_mask_path(std::string &fn, int k) {
  return fn + ".k" + std::to_string(k) + ".mask";
}

int ms_index(int argc, char *argv[]) {
  int c;
  bool usage = 0;
  int k = 0;
  std::vector<int> ls;
  while ((c = getopt(argc, argv, "k:l:h")) >= 0) {
    switch (c) {
    case 'k':
      k = atoi(optarg);
      break;
    case 'l':
      ls.push_back(atoi(optarg));
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
  auto ms = read_masked_superstring(fn);
  // If k is not set, infer it assuming the standard format of the mask.
  if (!k)
    k = infer_k(ms.mask);
  write_superstring(superstring_path, ms.superstring);
  // Construct and dump the BW-transformed masks.
  auto bw_transformed_masks = construct_bw_transformed_masks(ms, k, ls);
  for (auto [m, l] : bw_transformed_masks) {
    mask_dump(compute_mask_path(fn, l), m);
  }
  // Construct and dump the FM-index.
  fm_index_t fm_index;
  sdsl::construct(fm_index, superstring_path, 1);
  sdsl::store_to_file(fm_index, fn + ".fm9");
  // Clean the not needed superstring file.
  std::filesystem::remove(superstring_path);
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
  std::string kmer = argv[optind + 1];
  std::string index_path = fn + ".fm9";
  std::string mask_path = compute_mask_path(fn, (int)kmer.size());
  if (!std::filesystem::exists(std::filesystem::path{mask_path}) || !std::filesystem::exists(std::filesystem::path{index_path})) {
      std::cerr << "The index for file " << fn << " and k=" << std::to_string(kmer.size()) << " is not properly created." << std::endl;
      std::cerr << "Please run `./ms-index index` before." << std::endl;
      usage_query();
      return 1;
  }
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

int ms_clean(int argc, char *argv[]) {
    if (argc != 2 || strcmp(argv[1], "-h") == 0) {
        return usage_clean();
    }
    std::string fn = argv[1];
    std::filesystem::remove(fn + ".sstr");
    std::filesystem::remove(fn + ".fm9");
    for (int k = 1; k < 64; ++k) {
        std::filesystem::remove(compute_mask_path(fn, k));
    }
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
  else if (strcmp(argv[1], "clean") == 0)
      ret = ms_clean(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-v") == 0)
    return version();
  else
    return usage();

  return ret;
}
