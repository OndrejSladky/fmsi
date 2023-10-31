#include <filesystem>
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
  std::cerr
      << "MS-Index is a tool for efficient indexing of Masked Superstrings."
      << std::endl;
  std::cerr << std::endl << "The recognized commands are:" << std::endl;
  std::cerr << "  `index` - Creates a BWT based index of the given masked "
               "superstring."
            << std::endl;
  std::cerr << "  `query` - Queries a $k$-mer against an index." << std::endl;
  std::cerr << "  `clean` - Cleans the files stored for index." << std::endl;
  std::cerr << "  `-v`    - Prints the version of the program." << std::endl;
  std::cerr << "  `-h`    - Prints this help." << std::endl;
  return 1;
}

static int usage_index() {
  std::cerr
      << "MS-Index Index creates the index for a given masked superstring."
      << std::endl;
  std::cerr << std::endl << "The recognized arguments are:" << std::endl;
  std::cerr << "  `-p path_to_fasta` - The path to the fasta file with masked "
               "superstring to be indexed. This is a required argument."
            << std::endl;
  std::cerr << "  `-k value_of_k`    - The size of one k-mer. If not provided "
               "k is computed from the masked superstring under the assumption "
               "that the last run of zeros has length k-1."
            << std::endl;
  std::cerr << "  `-l value_of_l`    - The size of l for which a mask is "
               "computed from the mask for k. l should not be greater than k. "
               "This argument can be provided multiple times."
            << std::endl;
  std::cerr << "  `-h`               - Prints this help and terminates."
            << std::endl;
  return 1;
}

static int usage_query() {
  std::cerr << "MS-Index Query return whether the provided $k$-mer is in the "
               "masked superstring or not."
            << std::endl;
  std::cerr
      << "`./ms-index index` must be run on the provided fasta file beforehand."
      << std::endl;
  std::cerr << std::endl << "The recognized arguments are:" << std::endl;
  std::cerr << "  `-p path_to_fasta` - The path to the fasta file from which "
               "the index was created."
            << std::endl;
  std::cerr << "  `-f function`      - A function to determine whether a "
               "$k$-mer is represented based on the number of set and unset "
               "occurrences. The recognized functions are following:"
            << std::endl;
  std::cerr << "      `or`  - Consider $k$-mer represented when any of its "
               "occurrences is set."
            << std::endl;
  std::cerr << "      `all` - Assume that all occurrence are either set or "
               "unset and determine the presence by arbitrary occurrence."
            << std::endl;
  std::cerr << "      `and` - Consider $k$-mer represented when all its "
               "occurrences are set."
            << std::endl;
  std::cerr << "      `xor` - Consider $k$-mer represented when an odd number "
               "of occurrences is set."
            << std::endl;
  std::cerr << "      `X-Y` (where X and Y can be any integers) - Consider "
               "$k$-mer represented when its number of set occurrences is "
               "between X and Y (inclusive)."
            << std::endl;
  std::cerr << "  `-h` - Prints this help and terminates." << std::endl;
  std::cerr << std::endl << "The last positional argument is the queried k-mer" << std::endl;
  return 1;
}

static int usage_clean() {
  std::cerr
      << "MS-Index Index creates the index for a given masked superstring."
      << std::endl;
  std::cerr << std::endl << "The recognized arguments are:" << std::endl;
  std::cerr << "  `-p path_to_fasta` - The path to the fasta file with masked "
               "superstring that is indexed. This is a required argument."
            << std::endl;
  std::cerr << "  `-h`               - Prints this help and terminates."
            << std::endl;
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
  std::string fn;
  std::vector<int> ls;
  while ((c = getopt(argc, argv, "p:k:l:h")) >= 0) {
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
    case 'p':
      fn = optarg;
      break;
    default:
      return usage_index();
    }
  }
  if (usage) {
    usage_index();
    return 0;
  }
  if (fn == "") {
    std::cerr << "Path to the fasta file is a required argument." << std::endl;
    return usage_index();
  }
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
  std::string fn;
  std::function<bool(size_t, size_t)> f = nullptr;
  while ((c = getopt(argc, argv, "p:f:h")) >= 0) {
    switch (c) {
    case 'f':
      f = mask_function(optarg);
      break;
    case 'h':
      usage = true;
      break;
    case 'p':
      fn = optarg;
      break;
    default:
      return usage_query();
    }
  }
  if (usage) {
    usage_query();
    return 0;
  }
  if (fn == "") {
    std::cerr << "Path to the fasta file is a required argument." << std::endl;
    return usage_query();
  }
  if (optind + 1 > argc) {
    usage_query();
    return 1;
  }
  std::string kmer = argv[optind];
  std::string index_path = fn + ".fm9";
  std::string mask_path = compute_mask_path(fn, (int)kmer.size());
  if (!std::filesystem::exists(std::filesystem::path{mask_path}) ||
      !std::filesystem::exists(std::filesystem::path{index_path})) {
    std::cerr << "The index for file " << fn
              << " and k=" << std::to_string(kmer.size())
              << " is not properly created." << std::endl;
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
  bool usage = false;
  int c;
  std::string fn;
  while ((c = getopt(argc, argv, "p:f:h")) >= 0) {
    switch (c) {
    case 'h':
      usage = true;
      break;
    case 'p':
      fn = optarg;
      break;
    default:
      return usage_clean();
    }
  }
  if (usage) {
    usage_clean();
    return 0;
  }
  if (fn == "") {
    std::cerr << "Path to the fasta file is a required argument." << std::endl;
    return usage_clean();
  }
  std::filesystem::remove(fn + ".sstr");
  std::filesystem::remove(fn + ".fm9");
  for (int k = 1; k < 64; ++k) {
    std::filesystem::remove(compute_mask_path(fn, k));
  }
  return 0;
}

int main(int argc, char *argv[]) {
  int ret;
  if (argc < 2)
    return usage();
  if (strcmp(argv[1], "index") == 0)
    ret = ms_index(argc - 1, argv + 1);
  else if (strcmp(argv[1], "query") == 0)
    ret = ms_query(argc - 1, argv + 1);
  else if (strcmp(argv[1], "clean") == 0)
    ret = ms_clean(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-v") == 0)
    return version();
  else if (strcmp(argv[1], "-h") == 0) {
    usage();
    return 0;
  } else
    return usage();

  return ret;
}
