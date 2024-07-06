#pragma once

#include <vector>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v5.hpp>
#include "QSufSort.h"
#include "functions.h"

typedef unsigned char byte;

struct fms_index {
    sdsl::bit_vector ac_gt;
    sdsl::rank_support_v5<1> ac_gt_rank;
    sdsl::bit_vector ac;
    sdsl::rank_support_v5<1> ac_rank;
    sdsl::bit_vector gt;
    sdsl::rank_support_v5<1> gt_rank;
    sdsl::rrr_vector<> sa_transformed_mask;
    std::vector<size_t> counts;
    size_t dollar_position;
};

inline byte char_to_int(char c) {
    switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'a':
            return 0;
        case 'c':
            return 1;
        case 'g':
            return 2;
        case 't':
            return 3;
        default:
            return -1;
    }
}

inline bool is_upper(char c) {
    return c >= 'A' && c <= 'Z';
}

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

bool query(const fms_index& index, const std::string& pattern, demasking_function_t f) {
    size_t i = 0, i_rev = 0;
    size_t j = index.sa_transformed_mask.size(), j_rev = index.sa_transformed_mask.size();
    for (size_t k = pattern.size(); k > 0; --k) {
        update_range(index, i, j, char_to_int(pattern[k-1]));
    }
    // Separately optimize or.
    if (f == nullptr) {
        for (size_t k = i; k < j; ++k) {
            if (index.sa_transformed_mask[k]) {
                return true;
            }
        }
    }
    for (char p : pattern) {
        update_range(index, i_rev, j_rev, 3 ^ char_to_int(p));
    }
    // Separately optimize or.
    if (f == nullptr) {
        for (size_t k = i_rev; k < j_rev; ++k) {
            if (index.sa_transformed_mask[k]) {
                return true;
            }
        }
        return false;
    }
    // Do not optimize code for k-mers that are their own reverse complement as they're not very common.
    bool own_rc = true;
    for (size_t k = 0; k < pattern.size(); ++k) {
        if (char_to_int(pattern[k]) != (3 ^ char_to_int(pattern[pattern.size() - k - 1]))) {
            own_rc = false;
            break;
        }
    }
    size_t ones = 0;
    for (size_t k = i; k < j; ++k) {
        ones += index.sa_transformed_mask[k];
    }
    if (!own_rc) for (size_t k = i_rev; k < j_rev; ++k) {
        ones += index.sa_transformed_mask[k];
    }
    size_t total = j - i + j_rev - i_rev;
    if (own_rc) total = j - i;
    return f((int)ones, (int)total);
}

qsint_t* convert_superstring(std::string ms) {
    auto ret = new qsint_t[ms.size() + 1];
    for (size_t i = 0; i < ms.size(); ++i) {
        ret[i] = char_to_int(ms[i]);
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
        bwt[isa[i+1]] = char_to_int(ms[i]);
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
    index.ac_gt_rank = sdsl::rank_support_v5<1>(&index.ac_gt);
    index.ac_rank = sdsl::rank_support_v5<1>(&index.ac);
    index.gt_rank = sdsl::rank_support_v5<1>(&index.gt);

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
    index.ac_gt_rank = sdsl::rank_support_v5<1>(&index.ac_gt);
    sdsl::load_from_file(index.ac, basename + ".ac");
    index.ac_rank = sdsl::rank_support_v5<1>(&index.ac);
    sdsl::load_from_file(index.gt, basename + ".gt");
    index.gt_rank = sdsl::rank_support_v5<1>(&index.gt);
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
