#pragma once

#include <vector>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/rank_support_v.hpp>
#include "QSufSort.h"
#include "functions.h"
#include "kmers.h"

typedef unsigned char byte;

struct fms_index {
    sdsl::bit_vector ac_gt;
    sdsl::rank_support_v<1> ac_gt_rank;
    sdsl::bit_vector ac;
    sdsl::rank_support_v<1> ac_rank;
    sdsl::bit_vector gt;
    sdsl::rank_support_v<1> gt_rank;
    sdsl::rrr_vector<> sa_transformed_mask;
    std::vector<size_t> counts;
    size_t dollar_position;
};



size_t rank(const fms_index& index, size_t i, byte c) {
    auto gt_position = index.ac_gt_rank(i);
    if (c >= 2) { // G/T
        auto t_position = index.gt_rank(gt_position);
        if (c == 2) { // G
            return gt_position - t_position;
        } else { // T
            return t_position;
        }
    } else { // A/C
        auto c_position = index.ac_rank(i - gt_position);
        if (c == 0) { // A
            // Ranks are offset by 1 compared to the indices.
            return i - gt_position - c_position - (i >= index.dollar_position + 1);
        } else { // C
            return c_position;
        }
    }
}

byte access(const fms_index& index, size_t i) {
    auto gt_position = index.ac_gt_rank(i);
    if (index.ac_gt[i]) {
        return 2 + index.gt[gt_position];
    } else {
        return index.ac[i - gt_position];
    }
}

void update_range(const fms_index& index, size_t& i, size_t& j, byte c) {
    if (j == i) return;
    auto count = index.counts[c];
    i = count + rank(index, i, c);
    j = count + rank(index, j, c);
}

template <typename T>
void update_range_with_pattern(const fms_index& index, size_t& sa_start, size_t& sa_end, T pattern, int k) {
    if (sa_start == size_t(-1)) {
        sa_start = 0;
        sa_end = index.sa_transformed_mask.size();
        // Find the SA coordinates of the forward pattern.
        for (int i = 0; i < k; ++i) {
            update_range(index, sa_start, sa_end, pattern & 3);
            pattern >>= 2;
        }
    } else {
        // Streaming query.
        // TODO: Update with klcp.
        throw std::invalid_argument("Streaming queries are not supported yet.");
        int shift = 2 * (k - 1);
        T mask = 0b11 << shift;
        update_range(index, sa_start, sa_end, (pattern & mask) >> shift);
    }
}

template <bool maximized_ones=false, typename T>
int single_query_or(const fms_index& index, T pattern, int k, size_t& sa_start, size_t& sa_end) {
    update_range_with_pattern(index, sa_start, sa_end, pattern, k);
    // Separately optimize all-or-nothing and or.
    if constexpr (maximized_ones) {
        if (sa_start != sa_end) return index.sa_transformed_mask[sa_start];
        return -1;
    } else {
        for (size_t i = sa_start; i < sa_end; ++i) {
            if (index.sa_transformed_mask[i]) {
                return true;
            }
        }
        return (sa_start == sa_end) ? -1 : false;
    }
}

template <bool maximized_ones=false, typename T>
int single_query_or(const fms_index& index, T pattern, int k) {
    size_t sa_start = -1, sa_end = -1;
    return single_query_or<maximized_ones>(index, pattern, k, sa_start, sa_end);
}

template <typename T>
std::pair<size_t, size_t> single_query_general(const fms_index& index, T pattern, int k, size_t& sa_start, size_t& sa_end) {
    update_range_with_pattern(index, sa_start, sa_end, pattern, k);
    size_t ones = 0;
    for (size_t i = sa_start; i < sa_end; ++i) {
        ones += index.sa_transformed_mask[i];
    }
    return {ones, sa_end - sa_start};
}

template <typename T>
std::pair<size_t, size_t> single_query_general(const fms_index& index, T pattern, int k) {
    size_t sa_start = -1, sa_end = -1;
    return single_query_general(index, pattern, k, sa_start, sa_end);
}

template <bool maximized_ones = false, typename T>
std::vector<bool> query_kmers_streaming_forward_or(const fms_index& index, char* sequence, size_t sequence_length, int k) {
    T kmer = 0;
    for (int i = 0; i < k; ++i) {
        kmer = (kmer << 2) | nucleotideToInt[(uint8_t)sequence[i]];
    }
    std::vector<bool> result (sequence_length - k + 1);
    size_t sa_start = -1, sa_end = -1;
    result[0] = single_query_or<maximized_ones>(index, kmer, k, sa_start, sa_end);
    // Robust mask working even for k power of 2.
    T kmers_mask = (1 << (2 * k - 1));
    kmers_mask = kmers_mask | (kmers_mask - 1);

    for (size_t i = k; i < sequence_length; ++i) {
        kmer = (kmer << 2) | nucleotideToInt[(uint8_t)sequence[i]];
        kmer &= kmers_mask;
        result[i - k + 1] = single_query_or<maximized_ones>(index, kmer, k, sa_start, sa_end);
    }
    return result;
}


enum class query_mode {
    orr,
    all,
    general,
};

template <query_mode mode, typename T>
size_t query_kmers_single_bidir(const fms_index& index, char* sequence, size_t sequence_length, int k, demasking_function_t f) {
    T kmer = 0;
    for (int i = 0; i < k - 1; ++i) {
        kmer = (kmer << 2) | nucleotideToInt[(uint8_t)sequence[i]];
    }
    size_t result = 0;
    for (size_t i = k - 1; i < sequence_length; ++i) {
        kmer = (kmer << 2) | nucleotideToInt[(uint8_t)sequence[i]];
        if constexpr (mode != query_mode::general) {
            int got;
            if constexpr (mode == query_mode::orr) {
                got = single_query_or<false>(index, kmer, k);
            } else {
                got = single_query_or<true>(index, kmer, k);
            }
            if constexpr (mode == query_mode::orr) {
                if (got != 1) {
                    got = single_query_or<false>(index, ReverseComplement(kmer, k), k);
                }
            } else {
                if (got == -1) {
                    got = single_query_or<true>(index, ReverseComplement(kmer, k), k);
                }
            }
            result += got == 1;
        } else {
            auto [ones, total] = single_query_general(index, kmer, k);
            T rev_kmer = ReverseComplement(kmer, k);
            if (rev_kmer != kmer) {
                auto [ones_rev, total_rev] = single_query_general(index, rev_kmer, k);
                ones += ones_rev;
                total += total_rev;
            }
            if (f(ones, total)) {
                ++result;
            }
        }
    }
    return result;
}

template <query_mode mode, bool has_klcp, typename T>
size_t query_kmers(const fms_index& index, char* sequence, size_t sequence_length, int k, demasking_function_t f = nullptr) {
    if constexpr (has_klcp && mode != query_mode::general) {
        std::vector<bool> result = query_kmers_streaming_forward_or<T>(index, sequence, sequence_length, k);
        ReverseComplementString(sequence, sequence_length);
        std::vector<bool> result_rev = query_kmers_streaming_forward_or<T>(index, sequence, sequence_length, k);
        size_t ret = 0;
        for (size_t i = 0; i < result.size(); ++i) {
            if (result[i] || result_rev[result.size() - i - 1]) {
                ++ret;
            }
        }
        return ret;
    } else {
        return query_kmers_single_bidir<mode, T>(index, sequence, sequence_length, k, f);
    }
}







qsint_t* convert_superstring(std::string ms) {
    auto ret = new qsint_t[ms.size() + 1];
    for (size_t i = 0; i < ms.size(); ++i) {
        ret[i] = nucleotideToInt[(uint8_t)ms[i]];
    }
    return ret;
}

fms_index construct(std::string ms) {
    qsint_t *isa = convert_superstring(ms);
    // TODO: find out the required size of workspace.
    auto workspace = new qsint_t[ms.size() + 1];
    QSufSortSuffixSort(isa, workspace, (qsint_t)ms.size(),3, 0, 0);
    delete[] workspace;

    fms_index index;
    sdsl::bit_vector sa_transformed_mask(ms.size() + 1);
    std::vector<byte> bwt (ms.size() + 1);
    for (size_t i = 0; i < ms.size(); ++i) {
        bwt[isa[i+1]] = nucleotideToInt[(uint8_t)ms[i]];
        sa_transformed_mask[isa[i]] = is_upper(ms[i]);
    }

    index.dollar_position = isa[0];
    delete[] isa;
    index.sa_transformed_mask = sdsl::rrr_vector<>(sa_transformed_mask);
    sa_transformed_mask.resize(0);

    index.ac_gt = sdsl::bit_vector(bwt.size());
    size_t gt_count = 0;
    for (size_t i = 0; i < bwt.size(); ++i) {
        bool is_gt = bwt[i] >= 2;
        gt_count += is_gt;
        index.ac_gt[i] = is_gt;
    }
    size_t ac_count = bwt.size() - gt_count;
    index.ac = sdsl::bit_vector(ac_count);
    index.gt = sdsl::bit_vector(gt_count);
    size_t a_count = 0;
    size_t g_count = 0;
    size_t ac_index = 0;
    size_t gt_index = 0;
    for (size_t i = 0; i < bwt.size(); ++i) {
        bool is_one = bwt[i] & 1;
        if (index.ac_gt[i] == 0) {
            index.ac[ac_index++] = is_one;
            a_count += !is_one;
        } else {
            index.gt[gt_index++] = is_one;
            g_count += !is_one;
        }
    }
    index.counts = {1, a_count, ac_count, ac_count + g_count};
    index.ac_gt_rank = sdsl::rank_support_v<1>(&index.ac_gt);
    index.ac_rank = sdsl::rank_support_v<1>(&index.ac);
    index.gt_rank = sdsl::rank_support_v<1>(&index.gt);

    return index;
}

std::string export_ms(const fms_index& index) {
    std::string masked_letters = "acgtACGT";
    std::vector<char> ret(index.sa_transformed_mask.size() - 1);

    for (size_t i = 0, bw_index = 0; i < index.sa_transformed_mask.size() - 1; ++i) {
        byte letter = access(index, bw_index);
        bw_index = index.counts[letter] + rank(index, bw_index, letter);
        ret[index.sa_transformed_mask.size() - 2 - i] = masked_letters[letter + (index.sa_transformed_mask[bw_index] << 2)];
    }

    return {ret.begin(), ret.end()};
}

fms_index merge(const fms_index& a, const fms_index& b) {
    return construct(export_ms(a) + export_ms(b));
}

void dump_index(const fms_index& index, const std::string &fn) {
    auto basename = fn + ".fmsi";
    sdsl::store_to_file(index.ac_gt, basename + ".ac_gt");
    sdsl::store_to_file(index.ac, basename + ".ac");
    sdsl::store_to_file(index.gt, basename + ".gt");
    sdsl::store_to_file(index.sa_transformed_mask, basename + ".mask");
    std::ofstream out(basename + ".misc");
    out << index.dollar_position << std::endl;
    for (auto c : index.counts) {
        out << c << std::endl;
    }
    out.close();
}

fms_index load_index(const std::string &fn) {
    fms_index index;
    auto basename = fn + ".fmsi";
    sdsl::load_from_file(index.ac_gt, basename + ".ac_gt");
    index.ac_gt_rank = sdsl::rank_support_v<1>(&index.ac_gt);
    sdsl::load_from_file(index.ac, basename + ".ac");
    index.ac_rank = sdsl::rank_support_v<1>(&index.ac);
    sdsl::load_from_file(index.gt, basename + ".gt");
    index.gt_rank = sdsl::rank_support_v<1>(&index.gt);
    sdsl::load_from_file(index.sa_transformed_mask, basename + ".mask");
    std::ifstream in(basename + ".misc");
    in >> index.dollar_position;
    for (size_t i = 0; i < 4; ++i) {
        size_t c;
        in >> c;
        index.counts.push_back(c);
    }
    in.close();
    return index;
}
