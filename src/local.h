#pragma once

#include "parser.h"
#include "mask.h"
#include "functions.h"
#include "index.h"
#include "query.h"
#include <sdsl/rank_support_v5.hpp>

constexpr auto letters = "ACGT";

/// Convert the encoded KMer representation to string.
std::string number_to_k_mer(int encoded, int length) {
    std::string ret(length, 'N');
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret[length - i -1] = letters[encoded & 3];
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}

/// Set the mask such that the k-mer corresponding to the given interval is masked as ghost.
void delete_k_mer(sdsl::bit_vector &mask, size_t start, size_t count, size_t start_rc, size_t count_rc, assignable_function_t f) {
    bool zeroing = false;
    for (size_t i = 0; i < count; ++i) {
        if (f(i, count + count_rc) == 0) zeroing = true;
        mask[start + i] = !zeroing;
    }
    for (size_t i = 0; i < count_rc; ++i) {
        if (f(i + count, count + count_rc) == 0) zeroing = true;
        mask[start_rc + i] = !zeroing;
    }
}

/// Determine whether the given k-mer (or its RC) is represented.
/// If the k-mer appears in the superstring, determine its presence from
/// the total number of occurrences, the number of set occurrences via
/// the provided function f.
/// Delete the k-mer if present.
bool query_f_and_delete(fm_index_t &fm_index, sdsl::bit_vector &mask, sdsl::rank_support_v5<> &rank,
             assignable_function_t f, std::string kmer) {
    // TODO: refactor so that the code is deduplicated.
    size_t total = 0;
    size_t ones = 0;
    auto rc = reverse_complement(kmer);
    fm_index_t::size_type from, to;
    fm_index_t::size_type from_rc, to_rc;
    total += sdsl::backward_search(fm_index, 0, fm_index.size() - 1, kmer.begin(),
                                   kmer.end(), from, to);
    if (total) {
        for (size_t i = from; i <= to; ++i) {
            if (mask[i]) ones++;
        }
    }
    auto count_rc = sdsl::backward_search(fm_index, 0, fm_index.size() - 1,
                                          rc.begin(), rc.end(), from_rc, to_rc);
    if (count_rc) {
        total += count_rc;
        for (size_t i = from_rc; i <= to_rc; ++i) {
            if (mask[i]) ones++;
        }
    }
    bool res = f(ones, total);
    if (res) {
        delete_k_mer(mask, from, to - from + 1, from_rc, to_rc - from_rc + 1, f);
    }
    return res;
}

/// Find the extension to the provided last k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::string extension(fm_index_t &fm_index, sdsl::bit_vector &mask, sdsl::rank_support_v5<> &rank, std::string &base, assignable_function_t f, int d, bool left) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (int ext = 0; ext < (1 << (d << 1)); ++ext) {
        auto decoded_ext = number_to_k_mer(ext, d);
        auto next = left ? decoded_ext + base : base + decoded_ext;
        if (query_f_and_delete(fm_index, mask, rank, f, next)) {
            return decoded_ext;
        }
    }
    return "";
}

masked_superstring_t merge_superstrings(masked_superstring_t a, masked_superstring_t b) {
    std::vector<bool> mask;
    std::string superstring = a.superstring + b.superstring;
    for (auto m : a.mask) mask.push_back(m);
    for (auto m : b.mask) mask.push_back(m);
    assert(mask.size() == superstring.size());
    return {mask, superstring};
}

/// Find the next generalized simplitig.
masked_superstring_t next_generalized_simplitig(fm_index_t &fm_index, sdsl::bit_vector &mask, sdsl::rank_support_v5<> &rank, assignable_function_t f, int d_max, std::string &begin) {
    std::string last = begin, first = begin;
    int k = begin.size();
    std::string simplitig_back = begin;
    std::vector<bool> mask_back = {1};
    std::string simplitig_front = "";
    std::vector<bool> mask_front = {};
    int d_l = 1, d_r = 1;
    while (d_l <= d_max || d_r <= d_max) {
        if (d_r <= d_l) {
            auto base = last.substr(d_r);
            auto ext = extension(fm_index, mask, rank, base, f, d_r, false);
            if (ext.empty()) {
                // No right extension found.
                ++d_r;
            } else {
                // Extend the generalized simplitig to the right.
                last = base + ext;
                simplitig_back += ext;
                assert(ext.size() == d_r);
                for (int i = 1; i < d_r; ++i) mask_back.push_back(0);
                mask_back.push_back(1);
                d_r = 1;
                assert(mask_back.size() + k - 1 == simplitig_back.size());
            }
        } else {
            auto base = first.substr( 0, k - d_l);
            auto ext = extension(fm_index, mask, rank, base, f, d_l, true);
            if (ext.empty()) {
                // No left extension found.
                ++d_l;
            } else {
                // Extend the simplitig to the left.
                first = ext + base;
                std::reverse(ext.begin(), ext.end());
                simplitig_front += ext;
                assert(ext.size() == d_l);
                for (int i = 1; i < d_l; ++i) mask_front.push_back(0);
                mask_front.push_back(1);
                assert(mask_front.size() == simplitig_front.size());
                d_l = 1;
            }
        }
    }
    for (int i = 1; i < k; ++i) mask_back.push_back(0);
    assert(mask_back.size() == simplitig_back.size());
    std::reverse(mask_front.begin(), mask_front.end());
    std::reverse(simplitig_front.begin(), simplitig_front.end());
    return merge_superstrings( {mask_front, simplitig_front}, {mask_back, simplitig_back});
}

std::string next_k_mer(fm_index_t &fm_index, sdsl::bit_vector &mask, sdsl::rank_support_v5<> &rank, klcp_t &klcp, assignable_function_t f, size_t &start, int k) {
    while(start < fm_index.size()) {
        size_t end = start;
        while (end + 1 < klcp.size() && klcp[end+1]) {
            end++;
        }
        end++;
        size_t index = fm_index[start];
        if (index + k - 1 >= fm_index.size()) {
            start = end;
            continue;
        }
        auto k_mer = sdsl::extract(fm_index, index, index + k - 1);
        start = end;
        if (query_f_and_delete(fm_index, mask, rank, f, k_mer)) {
            return k_mer;
        }
    }
    // No more k-mers.
    return "";
}


masked_superstring_t local(fm_index_t &fm_index, sdsl::bit_vector &mask, sdsl::rank_support_v5<> &rank, klcp_t &klcp, assignable_function_t f, int k, int d_max) {
    size_t start = 0;
    masked_superstring_t ret;
    while (start < fm_index.size()) {
        auto k_mer = next_k_mer(fm_index, mask, rank, klcp, f, start, k);
        if (k_mer.empty()) break;
        ret = merge_superstrings(ret, next_generalized_simplitig(fm_index, mask, rank, f, d_max, k_mer));
    }
    return ret;
}