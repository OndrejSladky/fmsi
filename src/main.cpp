#include <filesystem>
#include <sdsl/suffix_arrays.hpp>

#include "fms_index.h"
#include "parser.h"
#include "version.h"
#include "compact.h"

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int usage() {
  std::cerr << "FMSI is a tool for efficient indexing of Masked Superstrings."
            << std::endl;
  std::cerr << std::endl << "The recognized commands are:" << std::endl;
  std::cerr << "  `index` - Creates a BWT based index of the given masked "
               "superstring."
            << std::endl;
  std::cerr << "  `query` - Queries a $k$-mer against an index." << std::endl;
  std::cerr << "  `clean` - Cleans the files stored for index." << std::endl;
  std::cerr << "  `merge` - Merges several indices." << std::endl;
  std::cerr << "  `compact` - Compacts the given index." << std::endl;
  std::cerr << "  `export` - Export the underlying masked superstring." << std::endl;
  std::cerr << "  `-v`    - Prints the version of the program." << std::endl;
  std::cerr << "  `-h`    - Prints this help." << std::endl;
  return 1;
}

static int usage_index() {
  std::cerr
      << "FMSI Index creates the index for a given masked superstring."
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

static int usage_merge() {
  std::cerr << "FMSI Merge merges several indices." << std::endl;
  std::cerr << std::endl << "The recognized arguments are:" << std::endl;
  std::cerr << "  `-p path_to_fasta`  - The path to the fasta file which "
               "should be merged. Can be provided multiple times. It is "
               "expected that it appears at least twice."
            << std::endl;
  std::cerr << "  `-r path_of_result` - The path where the result should be "
               "stored. This is a required argument."
            << std::endl;
  std::cerr << "  `-k value_of_k`     - The size of one k-mer." << std::endl;
  std::cerr << "  `-h`                - Prints this help and terminates."
            << std::endl;
  return 1;
}

static int usage_query() {
  std::cerr << "FMSI Query return whether the provided $k$-mer is in the "
               "masked superstring or not."
            << std::endl;
  std::cerr
      << "`./fmsi index` must be run on the provided fasta file beforehand."
      << std::endl;
  std::cerr << std::endl << "The recognized arguments are:" << std::endl;
  std::cerr << "  `-p path_to_fasta` - The path to the fasta file from which "
               "the index was created. Required."
            << std::endl;
  std::cerr << "  `-q path_to_queries` - The path to the file with k-mers to "
               "query. Required."
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

static int usage_normalize() {
  std::cerr << "FMSI Compact compacts the given index so that it "
               "does not occupy more space than needed."
            << std::endl;
  std::cerr << std::endl << "The recognized arguments are:" << std::endl;
  std::cerr << "  `-p path_to_fasta` - The path to the fasta file from which "
               "the index was created. Required."
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
  std::cerr << "  `-l` - Use local algorithm instead of global." << std::endl;
  std::cerr << "  `-d` - Value of d_max. Default 5." << std::endl;
  std::cerr << "  `-s` - Only print the masked superstring and do not "
               "compact the index"
            << std::endl;
  std::cerr << "  `-h` - Prints this help and terminates." << std::endl;
  return 1;
}

static int usage_export() {
    std::cerr << "FMSI Export exports the indexed masked superstring to its string form."
              << std::endl;
    std::cerr << std::endl << "The recognized arguments are:" << std::endl;
    std::cerr << "  `-p path_to_fasta` - The path to the fasta file from which "
                 "the index was created. Required."
              << std::endl;
    std::cerr << "  `-h` - Prints this help and terminates." << std::endl;
    return 1;
}

static int usage_clean() {
  std::cerr
      << "FMSI Index creates the index for a given masked superstring."
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

  std::cerr << "Starting " << fn << std::endl;
  std::string superstring_path = fn + ".sstr";
  auto ms = read_masked_superstring(fn);
  if (ms.size() == 0) {
    std::cerr << "The file '" << fn
              << "' is in incorrect format. It is supposed to be a fasta file "
                 "with a single entry, the masked superstring"
              << std::endl;
    return usage_index();
  }
  std::cerr << "Read masked superstring" << std::endl;
  // TODO: If k is not set, infer it assuming the standard format of the mask.
  fms_index index = construct(ms);
  std::cerr << "Constructed index" << std::endl;
  dump_index(index, fn);
  std::cerr << "Written index" << std::endl;
  return 0;
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
    std::cerr << "Path to the file with queries is a required argument."
              << std::endl;
    return usage_query();
  }

  fms_index index = load_index(fn);

  std::ifstream query_file(query_fn);
  std::string kmer;
  while (query_file >> kmer) {
    if (kmer.size() != size_t(k)) {
      std::cerr << "Skipped - size of queried k-mer " << kmer.size()
                << " is not k=" << k << std::endl;
      std::cout << "SKIPPED" << std::endl;
      continue;
    }
    bool found = query(index, kmer, f);
    if (found)
      std::cout << "FOUND\n";
    else
      std::cout << "NOT FOUND\n";
  }
  return 0;
}

int ms_merge(int argc, char *argv[]) {
  bool usage = false;
  int c;
  std::vector<std::string> fns;
  std::string result_fn;
  while ((c = getopt(argc, argv, "p:hr:")) >= 0) {
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
    default:
      return usage_merge();
    }
  }
  if (usage) {
    usage_merge();
    return 0;
  }

  fms_index res = load_index(fns[0]);
  std::cerr << "Loaded index " << fns[0] << std::endl;

  for (size_t i = 1; i < fns.size(); ++i) {
      res = merge(res, load_index(fns[i]));
      std::cerr << "Loaded and merged index " << fns[i] << std::endl;
  }

  dump_index(res, result_fn);
  std::cerr << "Result written" << std::endl;
  return 0;
}

int ms_normalize(int argc, char *argv[]) {
  bool usage = false;
  int c;
  int k = 0;
  int d_max = 5;
  bool only_print = false;
  bool use_local = false;
  std::string fn;
  std::function<bool(size_t, size_t)> f = mask_function("or");
  while ((c = getopt(argc, argv, "p:hk:d:f:sl")) >= 0) {
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
    case 'k':
      k = atoi(optarg);
      break;
    case 'd':
      d_max = atoi(optarg);
      break;
    case 's':
      only_print = true;
      break;
    case 'l':
      use_local = true;
      break;
    default:
      return usage_normalize();
    }
  }
  if (usage) {
    usage_normalize();
    return 0;
  }

  std::cerr << "Starting " << fn << std::endl;
  fms_index index = load_index(fn);
  std::cerr << "Loaded index" << std::endl;
  auto ms = export_ms(index);
  ms = normalize(ms, k, f);
  std::cerr << "Compacted" << std::endl;
  if (only_print) {
    std::cout << ">exported f-masked superstring" << std::endl;
    std::cout << ms << std::endl;
    return 0;
  }
  dump_index(construct(ms), fn);
  std::cerr << "Written index" << std::endl;
  return 0;
}

int ms_export(int argc, char *argv[]) {
    bool usage = false;
    int c;
    std::string fn;
    while ((c = getopt(argc, argv, "p:h")) >= 0) {
        switch (c) {
            case 'h':
                usage = true;
                break;
            case 'p':
                fn = optarg;
                break;
            default:
                return usage_export();
        }
    }
    if (usage) {
        usage_export();
        return 0;
    }

    fms_index index = load_index(fn);
    std::cout << ">exported f-masked superstring" << std::endl;
    std::cout << export_ms(index) << std::endl;
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

  std::filesystem::remove(fn + ".fmsi.ac_gt");
  std::filesystem::remove(fn + ".fmsi.ac");
  std::filesystem::remove(fn + ".fmsi.gt");
  std::filesystem::remove(fn + ".fmsi.mask");
  std::filesystem::remove(fn + ".fmsi.misc");
  std::cerr << "Cleaned " << fn << std::endl;
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
  // Recognize "normalize" for backwards compatibility
  else if (strcmp(argv[1], "normalize") == 0 || strcmp(argv[1], "compact") == 0)
    ret = ms_normalize(argc - 1, argv + 1);
  else if (strcmp(argv[1], "export") == 0)
      ret = ms_export(argc - 1, argv + 1);
  else if (strcmp(argv[1], "-v") == 0)
    return version();
  else if (strcmp(argv[1], "-h") == 0) {
    usage();
    return 0;
  } else
    return usage();

  return ret;
}
