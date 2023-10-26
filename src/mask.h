#pragma once

#include "QSufSort.h"
#include <fstream>
#include <iostream>
#include <sdsl/rrr_vector.hpp>
#include <vector>

typedef std::vector<bool> mask_t;
typedef sdsl::rrr_vector<63> bw_mask_t;
typedef sdsl::rank_support_rrr<1, 63> bw_mask_rank_t;

/// Compute the mask in suffix array coordinates.
bw_mask_t bw_transform_mask(const qsint_t *inverse_suffix_array,
                            mask_t original_mask) {
  sdsl::bit_vector result = sdsl::bit_vector(original_mask.size() + 1);
  for (uint64_t i = 0; i < original_mask.size(); ++i) {
    result[inverse_suffix_array[i]] = original_mask[i];
  }
  return bw_mask_t(result);
}

/// Serialize the mask to a given file.
void mask_dump(std::string fn, bw_mask_t mask) {
  std::ofstream of;
  of.open(fn, std::ios::out);
  mask.serialize(of);
  of.close();
}

/// Construct the mask from a given file.
bw_mask_t mask_restore(std::string fn) {
  std::ifstream f;
  f.open(fn, std::ios::in);
  bw_mask_t ret;
  ret.load(f);
  f.close();
  return ret;
}
