#pragma once

#include "QSufSort.h"
#include <fstream>
#include <iostream>
#include <sdsl/rrr_vector.hpp>
#include <vector>

typedef std::vector<bool> mask_t;
typedef sdsl::rrr_vector<4> bw_mask_t;

/// Compute the mask in suffix array coordinates.
bw_mask_t bw_transform_mask(const qsint_t *inverse_suffix_array,
                            mask_t original_mask) {
  sdsl::bit_vector result = sdsl::bit_vector(original_mask.size() + 1);
  for (uint64_t i = 0; i < original_mask.size(); ++i) {
    result[inverse_suffix_array[i]] = original_mask[i];
  }
  return bw_mask_t(result);
}

void mask_dump(std::string fn, bw_mask_t mask) {
  // TODO make this more efficient.
  std::ofstream of;
  of.open(fn, std::ios::binary | std::ios::out);
  mask.serialize(of);
  of.close();
}

bw_mask_t mask_restore(std::string fn) {
  // TODO make this more efficient.
  std::ifstream f;
  f.open(fn, std::ios::in);
  bw_mask_t ret;
  ret.load(f);
  f.close();
  return ret;
}
