#pragma once

#include <vector>
#include <filesystem>
#include <sdsl/select_support_mcl.hpp>
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
    sdsl::bit_vector klcp;
    int k;
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

void extend_range_with_klcp(const fms_index& index, size_t& i, size_t& j) {
    while(index.klcp[j-1]) j++;
    while(index.klcp[i-1]) i--;

    //auto rank = index.klcp_rank(i);
    // Select is indexed from 1.
    //i = index.klcp_select(rank) + 1;
    //j = index.klcp_select(rank + 1) + 1;
}

void get_range_with_pattern(const fms_index& index, size_t &sa_start, size_t &sa_end, char* pattern, int k) {
    sa_start = 0;
    sa_end = index.sa_transformed_mask.size();
    // Find the SA coordinates of the forward pattern.
    for (int i = k-1; i >= 0; --i) {
        update_range(index, sa_start, sa_end, nucleotideToInt[(uint8_t)pattern[i]]);
    }
}

template <bool maximized_ones=false>
int infer_presence(const fms_index& index, size_t sa_start, size_t sa_end) {
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
        if (sa_start == sa_end) {
            return -1;
        } else {
            return false;
        }
    }
}

template <bool maximized_ones=false>
int single_query_or(const fms_index& index, char* pattern, int k) {
    size_t sa_start = -1, sa_end = -1;
    get_range_with_pattern(index, sa_start, sa_end, pattern, k);
    return infer_presence<maximized_ones>(index, sa_start, sa_end);
}

std::pair<size_t, size_t> single_query_general(const fms_index& index, char* pattern, int k) {
    size_t sa_start, sa_end;
    get_range_with_pattern(index, sa_start, sa_end, pattern, k);
    size_t ones = 0;
    for (size_t i = sa_start; i < sa_end; ++i) {
        ones += index.sa_transformed_mask[i];
    }
    return {ones, sa_end - sa_start};
}

template <bool maximized_ones = false>
size_t query_kmers_streaming(const fms_index& index, char* sequence, char* rc_sequence, size_t sequence_length, int k) {
    std::vector<signed char> result (sequence_length - k + 1);
    size_t sa_start = -1, sa_end = -1;
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        size_t i_back = sequence_length - k - i;
        if (sa_start == sa_end) {
            get_range_with_pattern(index, sa_start, sa_end, sequence + i_back, k);
        } else {
            extend_range_with_klcp(index, sa_start, sa_end);
            update_range(index, sa_start, sa_end, nucleotideToInt[(uint8_t)sequence[i_back]]);
        }
        result[i_back] = infer_presence<maximized_ones>(index, sa_start, sa_end);
    }
    sa_start = sa_end = -1;
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        if (result[i] == 1 || (result[i] == 0 && maximized_ones)) {
            // This position can be skipped for performance.
            sa_start = sa_end = -1;
            continue;
        }
        size_t i_back = sequence_length - k - i;
        if (sa_start == sa_end) {
            get_range_with_pattern(index, sa_start, sa_end, rc_sequence + i_back, k);
        } else {
            extend_range_with_klcp(index, sa_start, sa_end);
            update_range(index, sa_start, sa_end, nucleotideToInt[(uint8_t)rc_sequence[i_back]]);
        }
        signed char res = infer_presence<maximized_ones>(index, sa_start, sa_end);
        result[i] = std::max(result[i], res);
    }
    return std::count(result.begin(), result.end(), 1);
}

enum class query_mode {
    orr,
    all,
    general,
};

template <query_mode mode>
size_t query_kmers_single(const fms_index& index, char* sequence, char* rc_sequence, size_t sequence_length, int k, demasking_function_t f) {
    size_t result = 0;
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        char *kmer = sequence + i;
        char *rc_kmer = rc_sequence + (sequence_length - k - i);
        if constexpr (mode != query_mode::general) {
            int got;
            if constexpr (mode == query_mode::orr) {
                got = single_query_or<false>(index, kmer, k);
            } else {
                got = single_query_or<true>(index, kmer, k);
            }
            if constexpr (mode == query_mode::orr) {
                if (got != 1) {
                    got = single_query_or<false>(index,rc_kmer , k);
                }
            } else {
                if (got == -1) {
                    got = single_query_or<true>(index, rc_kmer , k);
                }
            }
            result += got == 1;
        } else {
            auto [ones, total] = single_query_general(index, kmer, k);
            if (!AreStringsEqual(kmer, rc_kmer, k)) {
                auto [ones_rev, total_rev] = single_query_general(index, rc_kmer, k);
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

template <query_mode mode>
size_t query_kmers(const fms_index& index, char* sequence, size_t sequence_length, int k, bool has_klcp, demasking_function_t f = nullptr) {
    char *rc_sequence = ReverseComplementString(sequence, sequence_length);
    size_t ret = 0;
    if (has_klcp && mode != query_mode::general) {
        ret = query_kmers_streaming<mode==query_mode::all>(index, sequence, rc_sequence, sequence_length, k);
    } else {
        ret = query_kmers_single<mode>(index, sequence, rc_sequence, sequence_length, k, f);
    }
    delete[] rc_sequence;
    return ret;
}


template <typename T>
sdsl::bit_vector construct_klcp(qsint_t *isa, std::string& ms, size_t k_minus_1) {
    std::vector<qsint_t> sa(ms.size() + 1);
    for (size_t i = 0; i < ms.size() + 1; ++i) {
        sa[isa[i]] = i;
    }
    std::vector<T> kmers (ms.size() - k_minus_1 + 1);
    T kmer = 0;
    T mask = (1 << (2 * k_minus_1 - 1));
    mask = mask | (mask - 1);
    for (size_t i = 0; i < k_minus_1 - 1; ++i) {
        kmer = (kmer << 2) | nucleotideToInt[(uint8_t)ms[i]];
    }
    for (size_t i = 0; i < ms.size() - k_minus_1 + 1; ++i) {
        kmer = (kmer << 2);
        kmer |= nucleotideToInt[(uint8_t)ms[i+k_minus_1-1]];
        kmer &= mask;
        kmers[i] = kmer;
    }
    sdsl::bit_vector klcp (ms.size() + 1, 0);
    for (size_t i = 0; i < ms.size(); ++i) {
        qsint_t sa_i = sa[i];
        qsint_t sa_i1 = sa[i+1];
        if ((size_t)sa_i > ms.size() - k_minus_1 || (size_t)sa_i1 > ms.size() - k_minus_1) {
            continue;
        }
        klcp[i] = kmers[sa_i] == kmers[sa_i1];
    }
    return klcp;
}




qsint_t* convert_superstring(std::string ms) {
    auto ret = new qsint_t[ms.size() + 1];
    for (size_t i = 0; i < ms.size(); ++i) {
        ret[i] = nucleotideToInt[(uint8_t)ms[i]];
    }
    return ret;
}

template <typename T>
fms_index construct(std::string ms, int k, bool use_klcp) {
    qsint_t *isa = convert_superstring(ms);
    // TODO: find out the required size of workspace.
    auto workspace = new qsint_t[ms.size() + 1];
    QSufSortSuffixSort(isa, workspace, (qsint_t)ms.size(),3, 0, 0);
    delete[] workspace;

    fms_index index;

    if (use_klcp) {
        index.klcp = construct_klcp<T>(isa, ms, k-1);
    }

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

    index.k = k;

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
    if (a.k <= 32)  {
        return construct<uint64_t>(export_ms(a) + export_ms(b), a.k, a.klcp.size() > 0);
    } else {
        return construct<__uint128_t>(export_ms(a) + export_ms(b), a.k, a.klcp.size() > 0);
    }
}

void dump_index(const fms_index& index, const std::string &fn) {
    auto basename = fn + ".fmsi";
    sdsl::store_to_file(index.ac_gt, basename + ".ac_gt");
    sdsl::store_to_file(index.ac, basename + ".ac");
    sdsl::store_to_file(index.gt, basename + ".gt");
    sdsl::store_to_file(index.sa_transformed_mask, basename + ".mask");
    sdsl::store_to_file(index.klcp, basename + ".klcp");
    std::ofstream out(basename + ".misc");
    out << index.dollar_position << std::endl;
    for (auto c : index.counts) {
        out << c << std::endl;
    }
    out << index.k << std::endl;
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
    if (std::filesystem::exists(basename + ".klcp")) {
        sdsl::load_from_file(index.klcp, basename + ".klcp");
    }
    std::ifstream in(basename + ".misc");
    in >> index.dollar_position;
    for (size_t i = 0; i < 4; ++i) {
        size_t c;
        in >> c;
        index.counts.push_back(c);
    }
    in >> index.k;
    in.close();
    return index;
}
