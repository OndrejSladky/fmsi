#pragma once

#include <string>
#include <unordered_set>
#include <vector>
#include <deque>
#include <cstdint>
#include <algorithm>
#include <fstream>

#include "kmers.h"
#include "khash_utils.h"


/// Find the right extension to the provided last k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
std::pair<kmer_t, kmer_t> RightExtension(kmer_t last, kh_S_t *kMers, kh_wrapper_t wrapper, int k, int d, bool complements) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (kmer_t ext = 0; ext < (kmer_t(1) << (d << 1)); ++ext) {
        kmer_t next = BitSuffix(last, k - d) << (d << 1) | ext;
        if (containsKMer(kMers, wrapper, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the left extension to the provided first k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
std::pair<kmer_t, kmer_t> LeftExtension(kmer_t first, kh_S_t *kMers, kh_wrapper_t wrapper, int k, int d, bool complements) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (kmer_t ext = 0; ext < (kmer_t(1) << (d << 1)); ++ext) {
        kmer_t next = ext << ((k - d) << 1) | BitPrefix(first, k, k - d);
        if (containsKMer(kMers, wrapper, next, k, complements)) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the next generalized simplitig.
/// Update the provided superstring and the mask.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers only contain one k-mer from a complementary pair.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void NextGeneralizedSimplitig(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t begin, std::ostream& of,  int k, int d_max, bool complements) {
     // Maintain the first and last k-mer in the simplitig.
    kmer_t last = begin, first = begin;
    std::list<char> simplitig {NucleotideAtIndex(first, k, 0)};
    eraseKMer(kMers, wrapper, last, k, complements);
    int d_l = 1, d_r = 1;
    while (d_l <= d_max || d_r <= d_max) {
        if (d_r <= d_l) {
            auto extension = RightExtension(last, kMers, wrapper, k, d_r, complements);
            kmer_t ext = extension.first;
            if (ext == kmer_t(-1)) {
                // No right extension found.
                ++d_r;
            } else {
                // Extend the generalized simplitig to the right.
                eraseKMer(kMers, wrapper, extension.second, k, complements);
                for (int i = 1; i < d_r; ++i) simplitig.emplace_back((char)std::tolower(NucleotideAtIndex(last, k, i)));
                simplitig.emplace_back(NucleotideAtIndex(last, k, d_r));
                last = extension.second;
                d_r = 1;
            }
        } else {
            auto extension = LeftExtension(first, kMers,  wrapper, k, d_l, complements);
            kmer_t ext = extension.first;
            if (ext == kmer_t(-1)) {
                // No left extension found.
                ++d_l;
            } else {
                // Extend the simplitig to the left.
                eraseKMer(kMers, wrapper, extension.second, k, complements);
                for (int i = d_l - 1; i > 0; --i) simplitig.emplace_front((char)std::tolower(NucleotideAtIndex(extension.second, k, i)));
                simplitig.emplace_front(NucleotideAtIndex(extension.second, k, 0));
                first = extension.second;
                d_l = 1;
            }
        }
    }
    for (int i = 1; i < k; ++i) simplitig.emplace_back((char)std::tolower(NucleotideAtIndex(last, k, i)));
    of << std::string(simplitig.begin(), simplitig.end());
    of.flush();
}

/// Get the approximated shortest superstring of the given k-mers using the local greedy algorithm.
///
/// This runs in O(n d_max ^ k), where n is the number of k-mers, but for practical uses it is faster than AC version.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
/// Warning: this will destroy kMers.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void Local(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t _, std::ostream& of, int k, int d_max, bool complements) {
    size_t lastIndex = 0;
    while(true) {
        kmer_t begin = nextKMer(kMers, _, lastIndex);
        // No more k-mers.
        if (begin == kmer_t(-1)) return;
        NextGeneralizedSimplitig(kMers, wrapper, begin, of,  k, d_max, complements);
    }
}
