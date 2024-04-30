#pragma once

#include "functions.h"
#include "fms_index.h"
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
void count_k_mers(camel::kh_S64_t *k_mers, std::string superstring, std::vector<bool> mask,
                  int k, demasking_function_t f) {
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
    camel::khiter_t it = camel::kh_get_OCC64(total_occ, canonical_k_mer);
    if (it == kh_end(total_occ)) {
        it = camel::kh_put_OCC64(total_occ, canonical_k_mer, &ret);
    } else {
        total = kh_val(total_occ, it);
    }
    kh_value(total_occ, it) = ++total;
      it = camel::kh_get_OCC64(ones_occ, canonical_k_mer);
    if (it == kh_end(ones_occ)) {
        it = camel::kh_put_OCC64(ones_occ, canonical_k_mer, &ret);
    } else {
        ones = kh_val(ones_occ, it);
    }
    ones += mask[i];
    kh_value(ones_occ, it) = ones;
    bool contained = containsKMer(k_mers, canonical_k_mer, k, true);
    if (f(ones, total) && !contained) {
      kh_put_S64(k_mers, canonical_k_mer, &ret);
    } else if (!f(ones, total) && contained) {
      eraseKMer(k_mers, canonical_k_mer, k, true);
    }
  }

  camel::kh_destroy_OCC64(total_occ);
  camel::kh_destroy_OCC64(ones_occ);
}

/// Return the masked superstring corresponding to the given masked-cased
/// representation.
std::pair<std::vector<bool>, std::string> separate_mask_and_superstring(std::string ms) {
  auto mask = std::vector<bool>(ms.size());
  auto superstring = std::string(ms.size(), 'N');
  for (size_t i = 0; i < ms.size(); ++i) {
    char c = ms[i];
    mask[i] = is_upper(c);
    superstring[i] = _toupper(c);
  }
  return {mask, superstring};
}

/// Greedily compute a masked superstring with the same represented set as the
/// input.
std::string normalize(std::string ms, int k, demasking_function_t f) {
  camel::kh_S64_t *k_mers = camel::kh_init_S64();
  auto [mask, superstring] = separate_mask_and_superstring(ms);
  count_k_mers(k_mers, superstring, mask, k, f);
  std::stringstream ss;
  auto k_mer_vec = kMersToVec(k_mers);

  camel::PartialPreSort(k_mer_vec, k);
  camel::Global(k_mer_vec, ss, k, true);

  kh_destroy_S64(k_mers);
  return ss.str();
}