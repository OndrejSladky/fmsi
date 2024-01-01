#pragma once

#include "QSufSort.h"
#include "query.h"
#include <fstream>
#include <iostream>
#include <sdsl/rrr_vector.hpp>
#include <vector>

typedef std::vector<bool> mask_t;
typedef sdsl::rrr_vector<63> bw_mask_t;
typedef sdsl::rank_support_rrr<1, 63> bw_mask_rank_t;
// TODO: find out what the constants mean.
typedef sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127>>, 512, 1024>
        fm_index_t;

/// Compute the mask in suffix array coordinates.
bw_mask_t bw_transform_mask(const qsint_t *inverse_suffix_array,
                            mask_t original_mask) {
  sdsl::bit_vector result = sdsl::bit_vector(original_mask.size() + 1);
  for (uint64_t i = 0; i < original_mask.size(); ++i) {
    result[inverse_suffix_array[i]] = original_mask[i];
  }
  return {result};
}

/// Restore the original mask.
bw_mask_t bw_transform_mask(fm_index_t fm_index, mask_t mask) {
    auto bw_mask = sdsl::bit_vector(mask.size() + 1);
    for (size_t i = 0; i < mask.size(); ++i) {
        size_t index = fm_index[i];
        int value = index == mask.size() ? 0 : mask[index];
        bw_mask[i] = value;
    }
    return {bw_mask};
}

/// Restore the original mask.
mask_t bw_inverse_mask(fm_index_t fm_index, bw_mask_t bw_mask) {
    mask_t mask = mask_t(bw_mask.size() - 1);
    for (size_t i = 0; i < bw_mask.size(); ++i) {
        size_t index = fm_index[i];
        if (index < mask.size()) mask[fm_index[i]] = bw_mask[i];
    }
    return mask;
}

/// Copy elements from the source mask to the destination mask.
void merge_masks(mask_t &destination, mask_t source) {
    destination.reserve(destination.size() + source.size());
    for (auto x : source) destination.push_back(x);
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
