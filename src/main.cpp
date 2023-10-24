#include <sdsl/suffix_arrays.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "index.h"
#include "mask.h"
#include "parser.h"

// TODO: find out what the constants mean.
typedef sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024> fm_index_t;

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

int ms_index(int argc, char *argv[]) {
  int c;
  int usage = 0;
  int k = 13;
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
  write_superstring(superstring_path, ms.superstring);
  // Construct and dump the BW-transformed mask.
  mask_t bw_transformed_mask = construct_bw_transformed_mask(ms);
  mask_dump(mask_path, bw_transformed_mask);
  // Construct and dump the FM-index.
  fm_index_t fm_index;
  sdsl::construct(fm_index, superstring_path, 1);
  sdsl::store_to_file(fm_index, index_path);
  return 0;
}

int ms_query(int argc, char *argv[]) {
  if (optind + 1 > argc) {
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
  mask_t mask = mask_restore(mask_path);
  // Find the SA coordinates of the occurrences range of the k-mer and its RC.
  bool found = false;
  auto rc = reverse_complement(kmer);
  fm_index_t::size_type from, to;
  auto count = sdsl::backward_search(fm_index, 0, fm_index.size() - 1,
                                     kmer.begin(), kmer.end(), from, to);
  if (count) {
    found = mask[from];
  }
  auto count_rc = sdsl::backward_search(fm_index, 0, fm_index.size() - 1,
                                        rc.begin(), rc.end(), from, to);
  if (count_rc) {
    found |= mask[from];
  }

  if (found)
    printf("FOUND\n");
  else
    printf("NOT FOUND\n");
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
  else
    return usage();

  return ret;
}
