/*
  BW-transformed representation of mask.
  Author: Ondřej Sladký <ondra.sladky@gmail.com>
  Licence: MIT
*/

#pragma once

#include "QSufSort.h"
#include <fstream>
#include <iostream>
#include <vector>

typedef std::vector<bool> mask_t;

mask_t bw_transform_mask(const qsint_t *inverse_suffix_array,
                         mask_t original_mask) {
  mask_t result = mask_t(original_mask.size() + 1);
  for (uint64_t i = 0; i < original_mask.size(); ++i) {
    result[inverse_suffix_array[i]] = original_mask[i];
  }
  return result;
}

void mask_dump(const char *fn, mask_t mask) {
  // TODO make this more efficient.
  std::ofstream of;
  of.open(fn, std::ios::binary | std::ios::out);
  size_t size = mask.size();
  of.write((char *)&size, sizeof(size_t));
  for (size_t i = 0; i < size; ++i) {
    bool b = mask[i];
    of.write((char *)&b, sizeof(bool));
  }

  of.close();
}

mask_t mask_restore(const char *fn) {
  // TODO make this more efficient.
  std::ifstream f;
  f.open(fn, std::ios::binary | std::ios::in);
  size_t size;
  f.read((char *)&size, sizeof(size_t));
  auto ret = mask_t(size);
  for (size_t i = 0; i < size; ++i) {
    bool b;
    f.read((char *)&b, sizeof(bool));
    ret[i] = b;
  }

  f.close();
  return ret;
}
