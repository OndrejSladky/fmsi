#pragma once

#include "parser.h"
#include "QSufSort.h"
#include "mask.h"

#define ALPHABET_SIZE 4

inline qsint_t _convert_nucleotide(char c) {
  switch (c) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return 4;
  }
}

qsint_t *_convert_superstring(masked_superstring_t ms) {
  qsint_t *ret = (qsint_t*)malloc(ms.superstring.size());
  for (size_t i = 0; i < ms.superstring.size(); ++i) {
    ret[i] = _convert_nucleotide(ms.superstring[i]);
  }
  return ret;
}

mask_t construct_bw_transformed_mask(const char *fn, int k) {
  auto ms = read_masked_superstring(fn);
  append_reverse_complement(ms, k);
  qsint_t *superstring = _convert_superstring(ms);
  qsint_t *sa = (qsint_t*)malloc(ms.superstring.size());
  QSufSortSuffixSort(sa, superstring, (qsint_t)ms.superstring.size(),
                     (qsint_t)ALPHABET_SIZE - 1, 0, 0);
  return bw_transform_mask(sa, ms.mask);
}

#undef ALPHABET_SIZE