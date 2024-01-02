#pragma once

#include "QSufSort.h"
#include "compute_masks.h"
#include "mask.h"
#include "parser.h"
#include <sdsl/rrr_vector.hpp>
#include <sdsl/lcp_bitcompressed.hpp>

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

typedef std::vector<std::pair<bw_mask_t, int>> masks_with_k_t;

/// Return the mask indexed in the suffix array coordinates.
masks_with_k_t construct_bw_transformed_masks(masked_superstring_t ms, int k,
                                              std::vector<int> ls) {
  qsint_t *sa = _convert_superstring(ms);
  // TODO: find out the required size of workspace.
  qsint_t *workspace = new qsint_t[ms.superstring.size() + 1];
  QSufSortSuffixSort(sa, workspace, (qsint_t)ms.superstring.size(),
                     (qsint_t)ALPHABET_SIZE - 1, 0, 0);
  masks_with_k_t ret;
  ret.emplace_back(bw_transform_mask(sa, ms.mask), k);
  for (int l : ls)
    ret.emplace_back(bw_transform_mask(sa, compute_l_mask(ms.mask, k, l)), l);
  delete[] workspace;
  delete[] sa;
  return ret;
}

typedef sdsl::rrr_vector<63> klcp_t;

/// Construct the lcp array from the csa.
klcp_t construct_klcp(std::string fn, int k) {
    sdsl::lcp_bitcompressed<> lcp;
    sdsl::construct(lcp, fn, 1);
    sdsl::bit_vector klcp(lcp.size());
    for (size_t i = 0; i < lcp.size(); ++i) {
        klcp[i] = lcp[i] >= (unsigned)k;
    }
    return klcp_t(klcp);
}




/// Store the FM-index with its klcp array and masks to a associated file.
void dump_index_and_masks(std::string fn, fm_index_t fm_index,
                          masks_with_k_t bw_transformed_masks) {
  // Dump the masks.
  for (auto [m, l] : bw_transformed_masks) {
    // Include the k value in the name only if there are more masks.
    mask_dump(compute_mask_path(fn, l, bw_transformed_masks.size() > 1), m);
  }

  // Construct and dump the KLCP array.
  mask_dump(fn + ".klcp", construct_klcp(fn + ".sstr", bw_transformed_masks[0].second));

  // Dump the FM-index.
  sdsl::store_to_file(fm_index, fn + ".fm9");
}
