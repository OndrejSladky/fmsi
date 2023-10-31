#pragma once

#include "QSufSort.h"
#include "compute_masks.h"
#include "mask.h"
#include "parser.h"

constexpr int ALPHABET_SIZE = 4;

/// Return the integer representation of the nucleotide.
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

/// Return the integer representation of the given superstring to be used in
/// suffix sort.
qsint_t *_convert_superstring(masked_superstring_t ms) {
  size_t size = ms.superstring.size();
  qsint_t *ret = new qsint_t[size + 1];
  for (size_t i = 0; i < size; ++i) {
    ret[i] = _convert_nucleotide(ms.superstring[i]);
  }
  return ret;
}

/// Return the mask indexed in the suffix array coordinates.
std::vector<std::pair<bw_mask_t, int>>
construct_bw_transformed_masks(masked_superstring_t ms, int k,
                               std::vector<int> ls) {
  qsint_t *sa = _convert_superstring(ms);
  // TODO: find out the required size of workspace.
  qsint_t *workspace = new qsint_t[ms.superstring.size() + 1];
  QSufSortSuffixSort(sa, workspace, (qsint_t)ms.superstring.size(),
                     (qsint_t)ALPHABET_SIZE - 1, 0, 0);
  std::vector<std::pair<bw_mask_t, int>> ret;
  ret.emplace_back(bw_transform_mask(sa, ms.mask), k);
  for (int l : ls)
    ret.emplace_back(bw_transform_mask(sa, compute_l_mask(ms.mask, k, l)), l);
  delete[] workspace;
  delete[] sa;
  return ret;
}
