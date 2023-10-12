#pragma once

#include "QSufSort.h"
#include "mask.h"
#include "parser.h"

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
  size_t size = ms.superstring.size();
  // Allocate memory also for the
  qsint_t *ret = (qsint_t *)malloc((size + 1) * sizeof(qsint_t));
  for (size_t i = 0; i < size; ++i) {
    ret[i] = _convert_nucleotide(ms.superstring[i]);
  }
  return ret;
}

mask_t construct_bw_transformed_mask(const char *fn, int k) {
  auto ms = read_masked_superstring(fn);
  append_reverse_complement(ms, k);
  qsint_t *sa = _convert_superstring(ms);
  // TODO: find out the required size of workspace.
  qsint_t *workspace =
      (qsint_t *)malloc((ms.superstring.size() + 1) * sizeof(qsint_t));
  QSufSortSuffixSort(sa, workspace, (qsint_t)ms.superstring.size(),
                     (qsint_t)ALPHABET_SIZE - 1, 0, 0);
  auto ret = bw_transform_mask(sa, ms.mask);
  free(workspace);
  free(sa);
  return ret;
}

#undef ALPHABET_SIZE