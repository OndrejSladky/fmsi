#pragma once

#include "functions.h"
#include "mask.h"
#include "parser.h"

// Wrap up kmercamel in a namespace to avoid name conflicts.
namespace camel {
#include "kmercamel/src/global.h"
#include "kmercamel/src/khash_utils.h"
#include "kmercamel/src/kmers.h"
KHASH_MAP_INIT_INT64(OCC64, int)
} // namespace camel

#include <string>

/// Fill in kMers with the represented k-mers in the given superstring under f.
void count_k_mers(camel::kh_S64_t *k_mers, std::string superstring, mask_t mask,
                  int k, assignable_function_t f) {
  camel::kmer_t k_mer = 0;
  camel::kmer_t k_mer_mask = 1LL << (2 * k - 1);
  k_mer_mask |= k_mer_mask - 1;
  camel::kh_OCC64_t *total_occ = camel::kh_init_OCC64();
  camel::kh_OCC64_t *ones_occ = camel::kh_init_OCC64();

  for (int i = 0; i < k - 1; ++i) {
    k_mer = (k_mer << 2) | camel::NucleotideToInt(superstring[i]);
  }
  for (size_t i = 0; i < superstring.size() - k + 1; ++i) {
    k_mer = (k_mer << 2) | camel::NucleotideToInt(superstring[i + k - 1]);
    k_mer &= k_mer_mask;
    auto canonical_k_mer = std::min(k_mer, camel::ReverseComplement(k_mer, k));
    int total = 0, ones = 0;
    int ret;
    camel::khiter_t it = kh_put_OCC64(total_occ, canonical_k_mer, &ret);
    total = kh_val(total_occ, it);
    kh_value(total_occ, it) = total + 1;
    if (mask[i]) {
      it = kh_put_OCC64(ones_occ, canonical_k_mer, &ret);
      ones = kh_val(ones_occ, it);
      kh_value(ones_occ, it) = ones + 1;
    }
    bool contained = containsKMer(k_mers, k_mer, k, true);
    if (f(ones, total) && !contained) {
      kh_put_S64(k_mers, k_mer, &ret);
    } else if (!f(ones, total) && contained) {
      eraseKMer(k_mers, k_mer, k, true);
    }
  }

  camel::kh_destroy_OCC64(total_occ);
  camel::kh_destroy_OCC64(ones_occ);
}

/// Return the masked superstring corresponding to the given masked-cased representation.
masked_superstring_t separate_mask_and_superstring(std::string superstring) {
  masked_superstring_t ret;
  ret.mask = mask_t(superstring.size());
  ret.superstring = std::string(superstring.size(), 'N');
  for (size_t i = 0; i < superstring.size(); ++i) {
    char c = superstring[i];
    ret.mask[i] = _isupper(c);
    ret.superstring[i] = _toupper(c);
  }
  return ret;
}

/// Greedily compute a masked superstring with the same represented set as the input.
masked_superstring_t normalize(std::string superstring, mask_t mask, int k,
                               assignable_function_t f) {
  camel::kh_S64_t *k_mers = camel::kh_init_S64();
  count_k_mers(k_mers, superstring, mask, k, f);
  std::stringstream ss;
  auto k_mer_vec = kMersToVec(k_mers);

  camel::PartialPreSort(k_mer_vec, k);
  camel::Global(k_mer_vec, ss, k, true);

  kh_destroy_S64(k_mers);
  return separate_mask_and_superstring(ss.str());
}