#pragma once

#include "mask.h"

#include <stdexcept>

/// Compute an l mask from k mask that represents the set of l-mers
/// corresponding to the set of k-mer for the same superstring.
///
/// k must always be at least l.
mask_t compute_l_mask(mask_t &k_mask, int k, int l) {
  if (l > k)
    throw std::invalid_argument("l is larger than k when computing masks");
  mask_t ret(k_mask.size());
  int one_count = 0;
  for (size_t i = 0; i < k_mask.size(); ++i) {
    one_count += k_mask[i];
    int to_subtract = i - k + l - 1;
    if (to_subtract >= 0) {
      one_count -= k_mask[to_subtract];
    }
    ret[i] = bool(one_count);
  }
  return ret;
}

/// Infer k from the mask assuming the mask is in standard format -- the last
/// run of zeros has lenght k-1.
int infer_k(mask_t &mask) {
  int ret = 1;
  while (mask[mask.size() - size_t(ret)] == 0)
    ret++;
  return ret;
}