#pragma once

#include <functional>

#include "mask.h"
#include "parser.h"


/// Determine whether the given k-mer (or its RC) is represented.
/// If the k-mer appears in the superstring, determine its presence from
/// arbitrary occurrence.
bool query(fm_index_t fm_index, bw_mask_t mask, std::string kmer) {
  bool found = false;
  auto rc = reverse_complement(kmer);
  fm_index_t::size_type from, to;
  auto count = sdsl::backward_search(fm_index, 0, fm_index.size() - 1,
                                     kmer.begin(), kmer.end(), from, to);
  if (count) {
    // Determine the presence from the first occurrence in the SA order.
    found = mask[from];
  }
  auto count_rc = sdsl::backward_search(fm_index, 0, fm_index.size() - 1,
                                        rc.begin(), rc.end(), from, to);
  if (count_rc) {
    found |= mask[from];
  }
  return found;
}

/// Determine whether the given k-mer (or its RC) is represented.
/// If the k-mer appears in the superstring, determine its presence from
/// the total number of occurrences, the number of set occurrences via
/// the provided function f.
bool query_f(fm_index_t fm_index, bw_mask_t mask, bw_mask_rank_t rank,
             std::function<bool(int, int)> f, std::string kmer) {
  size_t total = 0;
  size_t ones = 0;
  auto rc = reverse_complement(kmer);
  fm_index_t::size_type from, to;
  total += sdsl::backward_search(fm_index, 0, fm_index.size() - 1, kmer.begin(),
                                 kmer.end(), from, to);
  if (total) {
    ones += rank(to + 1) - rank(from);
  }
  auto count_rc = sdsl::backward_search(fm_index, 0, fm_index.size() - 1,
                                        rc.begin(), rc.end(), from, to);
  if (count_rc) {
    total += count_rc;
    ones += rank(to + 1) - rank(from);
  }
  std::cerr << "Set occurrences: " << ones << std::endl
            << "Total occurrences: " << total << std::endl;
  return f(ones, total);
}