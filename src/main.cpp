#include <filesystem>
#include <sdsl/suffix_arrays.hpp>

#include "query.h"
#include "compute_masks.h"
#include "functions.h"
#include "index.h"
#include "version.h"
#include "parser.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fstream>

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
               "the index was created. Required."
            << std::endl;
    std::cerr << "  `-q path_to_queries` - The path to the file with k-mers to query. Required."
              << std::endl;
    std::cerr << "  `-k value_of_k` - The size of queried k-mers. Required."
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

int ms_index(int argc, char *argv[]) {
  int c;
  bool usage = false;
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
      usage = true;
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
  } else if (fn.empty()) {
    std::cerr << "Path to the fasta file is a required argument." << std::endl;
    return usage_index();
  }

  std::string superstring_path = fn + ".sstr";
  auto ms = read_masked_superstring(fn);
  // If k is not set, infer it assuming the standard format of the mask.
  if (!k)
    k = infer_k(ms.mask);
  write_superstring(superstring_path, ms.superstring);
  // Construct the BW-transformed masks.
  auto bw_transformed_masks = construct_bw_transformed_masks(ms, k, ls);
  // Construct the FM-index.
    fm_index_t fm_index;
    sdsl::construct(fm_index, superstring_path, 1);
    dump_index_and_masks(fn, fm_index, bw_transformed_masks);
    // Clean the not needed superstring file.
    std::filesystem::remove(superstring_path);
  return 0;
}

/// Restore the FM-index and the mask.
/// If the path is invalid, print error and return false. Otherwise return true and set the provided return arguments.
bool load_index_pair(std::string fn, int k, fm_index_t &ret_fm_index, bw_mask_t &ret_bw_mask) {
    std::string index_path = fn + ".fm9";
    std::string mask_path = compute_mask_path(fn, k, false);
    // If there is no general file for the mask, try k-specific.
    if (!std::filesystem::exists(std::filesystem::path{mask_path}))
        mask_path = compute_mask_path(fn, k, true);
    if (!std::filesystem::exists(std::filesystem::path{mask_path}) ||
        !std::filesystem::exists(std::filesystem::path{index_path})) {
        std::cerr << "The index for file " << fn
                  << " and k=" << std::to_string(k)
                  << " is not properly created." << std::endl;
        std::cerr << "Please run `./ms-index index` before." << std::endl;
        usage_query();
        return false;
    }
    // Load FM-index.
    sdsl::load_from_file(ret_fm_index, index_path);
    // Load mask.
    ret_bw_mask = mask_restore(mask_path);
    return true;
}

int ms_query(int argc, char *argv[]) {
    bool usage = false;
    int c;
    int k = 0;
    std::string fn;
    std::string query_fn;
    std::function<bool(size_t, size_t)> f = mask_function("or");
    while ((c = getopt(argc, argv, "p:f:hq:k:")) >= 0) {
        switch (c) {
            case 'f':
                try {
                    f = mask_function(optarg);
                } catch (std::invalid_argument &) {
                    std::cerr << "Function '" << optarg << "' not recognized." << std::endl;
                    return usage_query();
                }
                break;
            case 'h':
                usage = true;
                break;
            case 'p':
                fn = optarg;
                break;
            case 'q':
                query_fn = optarg;
                break;
            case 'k':
                k = atoi(optarg);
                break;
            default:
                return usage_query();
        }
    }
    if (usage) {
        usage_query();
        return 0;
    } else if (fn.empty()) {
        std::cerr << "Path to the fasta file is a required argument." << std::endl;
        return usage_query();
    } else if (query_fn.empty()) {
        std::cerr << "Path to the file with queries is a required argument." << std::endl;
        return usage_query();
    }


    fm_index_t fm_index;
    bw_mask_t mask;
    if (!load_index_pair(fn, k, fm_index, mask)) return 1;

    std::ifstream query_file(query_fn);
    std::string kmer;
    while(query_file >> kmer) {
        if (kmer.size() != size_t(k)) {
            std::cerr << "Skipped - size of queried k-mer " << kmer.size() << " is not k=" << k << std::endl;
            std::cout << "SKIPPED" << std::endl;
            continue;
        }
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
    }
    return 0;
}

int ms_merge(int argc, char *argv[]) {
  bool usage = false;
  int c;
  int k = 0;
  std::vector<std::string> fns;
  std::string result_fn;
  std::function<bool(size_t, size_t)> f = mask_function("or");
  while ((c = getopt(argc, argv, "p:hr:k:")) >= 0) {
    switch (c) {
    case 'h':
      usage = true;
      break;
    case 'p':
      fns.emplace_back(optarg);
      break;
    case 'r':
      result_fn = optarg;
      break;
        case 'k':
            k = atoi(optarg);
            break;
    default:
      return usage_query();
    }
  }
  if (usage) {
    usage_query();
    return 0;
  }

  std::string result_superstring;
  mask_t result_mask;
  for (auto& fn : fns) {
      fm_index_t fm_index;
      bw_mask_t mask;
      std::cerr << "Starting " << fn << std::endl;
      if (!load_index_pair(fn, k, fm_index, mask)) return 1;
      std::cerr << "Loaded " << fn << std::endl;

      merge_masks(result_mask, bw_inverse_mask(fm_index, mask));
      std::cerr << "Merged masks for " << fn << std::endl;
      result_superstring += sdsl::extract(fm_index, 0, fm_index.size() - 2);
      std::cerr << "Read superstring for " << fn << std::endl;
  }

  std::string result_superstring_fn = result_fn + ".sstr";
  write_superstring(result_superstring_fn, result_superstring);
  std::cerr << "Resulting superstring written" << std::endl;
    fm_index_t fm_index;
    sdsl::construct(fm_index, result_superstring_fn, 1);
    std::cerr << "Resulting FM-index constructed" << std::endl;

  masks_with_k_t masks = {{bw_transform_mask(fm_index, result_mask), k}};
    std::cerr << "Resulting mask constructed" << std::endl;

    dump_index_and_masks(result_fn, fm_index, masks);
    std::cerr << "Written mask and FM-index" << std::endl;
    // Clean the not needed superstring file.
    std::filesystem::remove(result_superstring_fn);
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
  std::filesystem::remove(compute_mask_path(fn, 0, false));
  for (int k = 1; k < 64; ++k) {
    std::filesystem::remove(compute_mask_path(fn, k, true));
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
  else if (strcmp(argv[1], "merge") == 0)
      ret = ms_merge(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-v") == 0)
    return version();
  else if (strcmp(argv[1], "-h") == 0) {
    usage();
    return 0;
  } else
    return usage();

  return ret;
}
