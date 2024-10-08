#pragma once

#include "global.h"
#include "kmers.h"

/// Return the length of the cycle cover which lower bounds the superstring length.
template <typename kmer_t, typename kh_wrapper_t>
size_t LowerBoundLength(kh_wrapper_t wrapper, std::vector<kmer_t> kMers, int k, bool complements) {
    auto cycle_cover = OverlapHamiltonianPath(wrapper, kMers, k, complements, true);
    size_t res = 0;
    for (auto &overlap : cycle_cover.second) {
        res += size_t(k) - size_t(overlap);
    }
    return res / (1 + complements);
}