/*
  BW-transformed representation of mask.
  Author: Ondřej Sladký <ondra.sladky@gmail.com>
  Licence: MIT
*/

#pragma once

#include "QSufSort.h"
#include <vector>

typedef std::vector<bool> mask_t;

mask_t bw_transform_mask(const qsint_t *suffix_array, mask_t original_mask) {
  mask_t result = mask_t(original_mask.size());
  for (uint64_t i = 0; i < original_mask.size(); ++i) {
    result[i] = original_mask[suffix_array[i]];
  }
  return result;
}

void mask_dump(const char *fn, const mask_t *mask);

void mask_restore(const char *fn, mask_t *mask);
