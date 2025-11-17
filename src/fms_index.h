#pragma once

#include <vector>
#include <filesystem>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include "QSufSort.h"
#include "functions.h"
#include "kmers.h"
#include <iostream>

typedef unsigned char byte;

/// A saturating counter that predicts whether the next queried k-mer is from the same strand as MS or from the reverse one.
struct strand_predictor {
    int score = 0;
    int* result_scores = new int[2] {0, 0};
    int previous = -1;

    inline int clipped(int x, int clipper = 7) {
        return std::max(-clipper, std::min(clipper, x));
    }

    void log_result (int forward_query_result, int reverse_query_result) {
        int difference = forward_query_result - reverse_query_result;
        score += difference;
        score = clipped(score);

        if (previous != -1) {
            result_scores[previous] += difference;
            result_scores[previous] = clipped(result_scores[previous]);
        }

        previous = forward_query_result > reverse_query_result;
    }

    /// Make prediction conditioned on the result of the previous query.
    bool predict_swap() {
        if (previous != -1 && result_scores[previous] != 0) {
            return result_scores[previous] < 0;
        } else {
            // If no information (score==0), assume forward strand.
            return score < 0;
        }
    }
};

constexpr int RRR_BLOCK_SIZE = 63;
struct fms_index {
    sdsl::bit_vector ac_gt;
    sdsl::rank_support_v5<1> ac_gt_rank;
    sdsl::bit_vector ac;
    sdsl::rank_support_v5<1> ac_rank;
    sdsl::bit_vector gt;
    sdsl::rank_support_v5<1> gt_rank;
    sdsl::rrr_vector<RRR_BLOCK_SIZE> sa_transformed_mask;
    sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE> mask_rank;
    std::vector<size_t> counts;
    size_t dollar_position;
    sdsl::bit_vector klcp;
    int k;
    strand_predictor predictor = strand_predictor();
    sdsl::select_support_hyb<1> ac_gt_select;
    sdsl::select_support_hyb<1> ac_select;
    sdsl::select_support_hyb<1> gt_select;
    sdsl::select_support_rrr<1> mask_select;
};


inline size_t rank(const fms_index& index, size_t i, byte c) {
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

inline byte access(const fms_index& index, size_t i) {
    auto gt_position = index.ac_gt_rank(i);
    if (index.ac_gt[i]) {
        return 2 + index.gt[gt_position];
    } else {
        return index.ac[i - gt_position];
    }
}

inline size_t select(const fms_index& index, const byte c, size_t i) {
    size_t pos;
    if (c >= 2) {
        pos = index.gt_select(i);
        if (c == 2) pos = i - pos;
    } else {
        pos = index.ac_select(i);
        if (c == 0) pos = i - pos;
    }
    size_t res = index.ac_gt_select(pos);
    if (c <= 1) res = pos - res;
    return res;
}

inline byte first_column_access(const fms_index& index, size_t i) {
    byte ret = -1;
    while (ret != 4 && i < index.counts[ret + 1]) {
        ret++;
    }
    return ret;
}

/// Go from range (i,j) for pattern P to range for c+P
inline void update_range(const fms_index& index, size_t& i, size_t& j, byte c) {
    if (j == i) return;
    auto count = index.counts[c];
    i = count + rank(index, i, c);
    j = count + rank(index, j, c);
}

/// Go from range (i,j) for pattern Px to range for P.
inline void extend_range_with_klcp(const fms_index& index, size_t& i, size_t& j) {
    while(index.klcp[j-1]) j++;
    while(index.klcp[i-1]) i--;

    // Unused alternative for klcp-based extending. Unused for being too memory heavy and slower.
    //auto rank = index.klcp_rank(i);
    // Select is indexed from 1.
    //i = index.klcp_select(rank) + 1;
    //j = index.klcp_select(rank + 1) + 1;
}

void get_range_with_pattern(const fms_index& index, size_t &sa_start, size_t &sa_end, char* pattern, int k) {
    sa_start = 0;
    sa_end = index.sa_transformed_mask.size();
    // Find the SA coordinates of the forward pattern.
    for (int i = k-1; i >= 0 && sa_start != sa_end; --i) {
        update_range(index, sa_start, sa_end, nucleotideToInt[(uint8_t)pattern[i]]);
    }
}

template <bool maximized_ones=false>
inline int infer_presence(const fms_index& index, size_t sa_start, size_t sa_end) {
    // Separately optimize all-or-nothing and or.
    if constexpr (maximized_ones) {
        if (sa_start != sa_end) return index.sa_transformed_mask[sa_start];
        return -1;
    } else {
        for (size_t i = sa_start; i < sa_end; ++i) {
            if (index.sa_transformed_mask[i]) {
                return 1;
            }
        }
        if (sa_start == sa_end) {
            return -1;
        } else {
            return 0;
        }
    }
}

inline int64_t kmer_order(const fms_index& index, size_t sa_start) {
    return index.mask_rank(sa_start);
}

inline int64_t kmer_order_if_present(const fms_index& index, size_t sa_start, size_t sa_end) {
    if (infer_presence<false>(index, sa_start, sa_end) == 1) {
        return kmer_order(index, sa_start);
    } else {
        return -1;
    }
}

template <bool maximized_ones=false>
int single_query_or(fms_index& index, char* pattern, int k) {
    size_t sa_start = -1, sa_end = -1;
    get_range_with_pattern(index, sa_start, sa_end, pattern, k);
    return infer_presence<maximized_ones>(index, sa_start, sa_end);
}

int64_t single_query_order(fms_index& index, char* pattern, int k) {
    size_t sa_start = -1, sa_end = -1;
    get_range_with_pattern(index, sa_start, sa_end, pattern, k);
    return kmer_order_if_present(index, sa_start, sa_end);
}

template <bool maximized_ones=false>
int64_t single_query_order_nonminimal(fms_index &index, char* pattern, int k) {
    size_t sa_start = -1, sa_end = -1;
    get_range_with_pattern(index, sa_start, sa_end, pattern, k);
    bool present = infer_presence<maximized_ones>(index, sa_start, sa_end);
    if (present >= 0) return sa_start;
    return -1;
}

std::pair<size_t, size_t> single_query_general(fms_index& index, char* pattern, int k) {
    size_t sa_start, sa_end;
    get_range_with_pattern(index, sa_start, sa_end, pattern, k);
    size_t ones = 0;
    for (size_t i = sa_start; i < sa_end; ++i) {
        ones += index.sa_transformed_mask[i];
    }
    return {ones, sa_end - sa_start};
}

template <bool maximized_ones = false>
void query_kmers_streaming(fms_index& index, char* sequence, char* rc_sequence, size_t sequence_length, int k, bool output_orders, std::ostream& of) {
    std::vector<int64_t> result (sequence_length - k + 1, -1);
    // Use saturating counter to ensure that RC strings are visited as forward strings.
    bool should_swap = index.predictor.predict_swap();
    if (should_swap) {
        std::swap(sequence, rc_sequence);
    }
    // Search on the forward strand.
    int forward_predictor_result = 0, backward_predictor_result = 0;
    size_t sa_start = -1, sa_end = -1;
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        size_t i_back = sequence_length - k - i;
        if (sa_start == sa_end) {
            get_range_with_pattern(index, sa_start, sa_end, sequence + i_back, k);
        } else {
            extend_range_with_klcp(index, sa_start, sa_end);
            update_range(index, sa_start, sa_end, nucleotideToInt[(uint8_t)sequence[i_back]]);
        }
        if (output_orders) {
            result[i_back] = kmer_order_if_present(index, sa_start, sa_end);
            if (result[i_back] >= 0) forward_predictor_result++;
            else forward_predictor_result--;
        } else {
            result[i_back] = infer_presence<maximized_ones>(index, sa_start, sa_end);
            forward_predictor_result += result[i_back];
        }
    }
    // Search on the reverse strand.
    sa_start = sa_end = -1;
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        if ((result[i] >= 0 && output_orders) || result[i] == 1 || (result[i] == 0 && maximized_ones)) {
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
        int64_t res;
        if (output_orders) {
            res = kmer_order_if_present(index, sa_start, sa_end);
            if (res >= 0) backward_predictor_result ++;
            else backward_predictor_result--;
        } else {
            res = infer_presence<maximized_ones>(index, sa_start, sa_end);
            backward_predictor_result += res;
        }
        result[i] = std::max(result[i], res);
    }
    // Log the results to the saturating counter for better future performance.
    if (should_swap) {
        std::reverse(result.begin(), result.end());
        std::swap(forward_predictor_result, backward_predictor_result);
    }
    index.predictor.log_result(forward_predictor_result, backward_predictor_result);
    
    for (size_t i = 0; i < result.size(); ++i) {
        int64_t c = result[i];
        if (output_orders) {
            if (i > 0) of << ",";
            of << c;
        }
        else if (c==1) {
            of << "1";
        } else {
            of << "0";
        }
    }
}


template <bool minimal_hash, bool maximized_ones>
void query_kmers_streaming_list(fms_index& index, char* sequence, char* rc_sequence, size_t sequence_length, int k, bool output_orders, std::vector<int64_t> result, size_t result_offset) {
    // Use saturating counter to ensure that RC strings are visited as forward strings.
    bool should_swap = index.predictor.predict_swap();
    if (should_swap) {
        std::swap(sequence, rc_sequence);
    }
    // Search on the forward strand.
    int forward_predictor_result = 0, backward_predictor_result = 0;
    size_t sa_start = -1, sa_end = -1;
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        size_t i_back = sequence_length - k - i;
        if (sa_start == sa_end) {
            get_range_with_pattern(index, sa_start, sa_end, sequence + i_back, k);
        } else {
            extend_range_with_klcp(index, sa_start, sa_end);
            update_range(index, sa_start, sa_end, nucleotideToInt[(uint8_t)sequence[i_back]]);
        }
        if (output_orders) {
            int64_t res;
            constexpr if (minimal_hash) {
                res = kmer_order_if_present(index, sa_start, sa_end);
            } else {
                res = sa_start
            }
            result[result_offset + i_back] = res;
            if (result[result_offset + i_back] >= 0) forward_predictor_result++;
            else forward_predictor_result--;
        } else {
            result[result_offset + i_back] = infer_presence<maximized_ones>(index, sa_start, sa_end);
            forward_predictor_result += result[result_offset + i_back];
        }
    }
    // Search on the reverse strand.
    sa_start = sa_end = -1;
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        if ((result[result_offset + i] >= 0 && output_orders) || result[result_offset + i] == 1 || (result[result_offset + i] == 0 && maximized_ones)) {
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
        int64_t res;
        if (output_orders) {
            constexpr if (minimal_hash) {
                res = kmer_order_if_present(index, sa_start, sa_end);
            } else {
                res = sa_start
            }
            if (res >= 0) backward_predictor_result++;
            else backward_predictor_result--;
        } else {
            res = infer_presence<maximized_ones>(index, sa_start, sa_end);
            backward_predictor_result += res;
        }
        result[i] = std::max(result[i], res);
    }
    // Log the results to the saturating counter for better future performance.
    if (should_swap) {
        std::reverse(result.begin(), result.end());
        std::swap(forward_predictor_result, backward_predictor_result);
    }
    index.predictor.log_result(forward_predictor_result, backward_predictor_result);
}


/// A copy of what happens in main. Unless chunking is done properly, this should be the main entry point. 
template <bool minimal_hash, bool maximized_ones>
std::vector<int64_t> query_kmers_streaming_with_chunking(fms_index& index, char* sequence, char* rc_sequence, size_t sequence_length, int k, bool output_orders) {
    std::vector<int64_t> result (sequence_length - k + 1, -1);
    size_t result_offset = 0;
    
    // Small overhead for the chunking (while gaining superior time from prediction).
    int64_t max_sequence_chunk_length = 400;
    max_sequence_chunk_length = k + std::max((int64_t)10, std::min(max_sequence_chunk_length, 2*(int64_t)std::sqrt(sequence_length)));

    while (sequence_length > 0) {
        int64_t current_length = next_invalid_character_or_end(sequence, sequence_length);
        
        while (current_length >= k) {
            int64_t chunk_length = std::min(current_length, max_sequence_chunk_length);
            query_kmers_streaming<minimal_hash, maximized_ones>(index, sequence, rc_sequence, sequence_length, k, output_orders, of);
            sequence += chunk_length - k + 1;
            current_length -= chunk_length - k + 1;
            sequence_length -= chunk_length - k + 1;
            result_offset += chunk_length - k + 1; 
        }
        // Skip also the next character.
        sequence_length -= current_length + 1;
        sequence += current_length + 1;
        // Print 0 on invalid k-mers.
        if (sequence_length >= 0) {
          for (int64_t i = 0; i < std::min((int64_t) k, current_length + 1); ++i) {
            if (output_orders) {
              result[result_offset++] = -1;
            }
            else {
              result[result_offset++] = 0;
            }
          }
        }
    }
    return result;
}

enum class query_mode {
    orr,
    all,
    general,
};

/// For each k-mer output 1 if it is found and 0 otherwise to the [of] stream.
template <query_mode mode>
void query_kmers_single(fms_index& index, char* sequence, char* rc_sequence, size_t sequence_length, int k, std::ostream& of, bool output_orders, demasking_function_t f) {
    for (size_t i = 0; i <= sequence_length - k; ++i) {
        char *kmer = sequence + i;
        char *rc_kmer = rc_sequence + (sequence_length - k - i);
        bool should_swap = index.predictor.predict_swap();
        if (should_swap) {
            std::swap(kmer, rc_kmer);
        }
        int forward_predictor_result = 0, backward_predictor_result = 0;
        if constexpr (mode != query_mode::general) {
            int64_t got;
            if (output_orders) {
                got = single_query_order(index, kmer, k);
            } else if constexpr (mode == query_mode::orr) {
                got = single_query_or<false>(index, kmer, k);
            } else {
                got = single_query_or<true>(index, kmer, k);
            }
            forward_predictor_result = got;
            if (output_orders) {
                if (forward_predictor_result >= 0) forward_predictor_result = 1;
                else {
                    got = single_query_order(index, rc_kmer, k);
                    backward_predictor_result = got >= 0 ? 1 : -1;
                }
            } else if constexpr (mode == query_mode::orr) {
                if (got != 1) {
                    got = single_query_or<false>(index,rc_kmer , k);
                    backward_predictor_result = got;
                }
            } else {
                if (got == -1) {
                    got = single_query_or<true>(index, rc_kmer , k);
                    backward_predictor_result = got;
                }
            }
            
            if (output_orders) {
                if (i > 0) of << ",";
                of << got;
            }
            else if (got == 1) {
                of << "1";
            } else {
                of << "0";
            }

            // Update strand predictor.
            if (should_swap) {
                std::swap(forward_predictor_result, backward_predictor_result);
            }
            index.predictor.log_result(forward_predictor_result, backward_predictor_result);
        } else {
            auto [ones, total] = single_query_general(index, kmer, k);
            // Do not count self complementary k-mers twice.
            if (!AreStringsEqual(kmer, rc_kmer, k)) {
                auto [ones_rev, total_rev] = single_query_general(index, rc_kmer, k);
                ones += ones_rev;
                total += total_rev;
            }
            if (f(ones, total)) {
                of << "1";
            } else {
                of << "0";
            }
        }
    }
}

template <query_mode mode>
void query_kmers(fms_index& index, char* sequence, size_t sequence_length, int k, bool has_klcp, std::ostream& of, bool output_orders, demasking_function_t f = nullptr) {
    char *rc_sequence = ReverseComplementString(sequence, sequence_length);
    if (has_klcp && mode != query_mode::general) {
        query_kmers_streaming<mode==query_mode::all>(index, sequence, rc_sequence, sequence_length, k, output_orders, of);
    } else {
        query_kmers_single<mode>(index, sequence, rc_sequence, sequence_length, k, of, output_orders, f);
    }
    delete[] rc_sequence;
}


template <typename T>
inline T obtain_kmer(std::vector<T> &kmers, std::string &ms, size_t i, int kmer_sparsity, T mask, int k_minus_1) {
    size_t i_base = i - (i % kmer_sparsity);
    T kmer = kmers[i_base / kmer_sparsity];
    for (size_t j = 0; j < i % kmer_sparsity; ++j) {
        kmer = (kmer << 2);
        kmer |= nucleotideToInt[(uint8_t)ms[i_base + k_minus_1 + j]];
        kmer &= mask;
    }
    return kmer;
}

template <typename T>
sdsl::bit_vector construct_klcp(qsint_t *sa, std::string& ms, size_t k_minus_1) {
    int kmer_sparsity = 4;
    std::vector<T> kmers ((ms.size() - k_minus_1 + 1) / kmer_sparsity + 1);
    T kmer = 0;
    T mask = (T(1) << (2 * k_minus_1 - 1));
    mask = mask | (mask - 1);
    for (size_t i = 0; i < k_minus_1 - 1; ++i) {
        kmer = (kmer << 2) | nucleotideToInt[(uint8_t)ms[i]];
    }
    for (size_t i = 0; i < ms.size() - k_minus_1 + 1; ++i) {
        kmer = (kmer << 2);
        kmer |= nucleotideToInt[(uint8_t)ms[i+k_minus_1-1]];
        kmer &= mask;
        if (i % kmer_sparsity == 0) {
            kmers[i / kmer_sparsity] = kmer;
        }
    }
    sdsl::bit_vector klcp (ms.size() + 1, 0);
    for (size_t i = 0; i < ms.size(); ++i) {
        qsint_t sa_i = sa[i];
        qsint_t sa_i1 = sa[i+1];
        if ((size_t)sa_i > ms.size() - k_minus_1 || (size_t)sa_i1 > ms.size() - k_minus_1) {
            continue;
        }
        klcp[i] = obtain_kmer(kmers, ms, sa_i, kmer_sparsity, mask, k_minus_1) == obtain_kmer(kmers, ms, sa_i1, kmer_sparsity, mask, k_minus_1);
    }
    return klcp;
}

auto _letters = {'A', 'C', 'G', 'T'};
template <bool minimal_hash>
char* kmer_access(const fms_index &index, int64_t identifier, int k) {
    char* result = (char*)(malloc(sizeof(char) * k));
    size_t position = identifier;
    constexpr if (minimal_hash) {
        position = index.mask_select(identifier);
    }

    for (int i = 0; i < k; ++i) {
        byte nucleotide = first_column_access(position);
        result[i] = _letters[nucleotide];
        size_t nucleotide_order = position - index.counts[nuclotide];
        position = select(index, nucleotide, nucleotide_order);
    }

    return result;
}


qsint_t* convert_superstring(std::string ms) {
    auto ret = new qsint_t[ms.size() + 1];
    for (size_t i = 0; i < ms.size(); ++i) {
        ret[i] = nucleotideToInt[(uint8_t)ms[i]];
    }
    return ret;
}

template <typename T>
fms_index construct(std::string &ms, int k, bool use_klcp) {
    qsint_t *qms = convert_superstring(ms);
    auto sa = new qsint_t[ms.size() + 1];
    QSufSortSuffixSort(qms, sa, (qsint_t)ms.size(),3, 0, 0);
    QSufSortGenerateSaFromInverse(qms, sa, (qsint_t)ms.size());
    delete[] qms;

    fms_index index;

    if (use_klcp) {
        index.klcp = construct_klcp<T>(sa, ms, k-1);
    }

    sdsl::bit_vector sa_transformed_mask(ms.size() + 1);
    std::vector<byte> bwt (ms.size() + 1);
    for (size_t i = 0; i <= ms.size(); ++i) {
        if (sa[i] == 0) {
            index.dollar_position = i;
        } else {
            bwt[i] = nucleotideToInt[(uint8_t)ms[sa[i] - 1]];
        }
        if (sa[i] != (qsint_t)ms.size()) {
            sa_transformed_mask[i] = is_upper(ms[sa[i]]);
        }
    }
    delete[] sa;
    index.sa_transformed_mask = sdsl::rrr_vector<RRR_BLOCK_SIZE>(sa_transformed_mask);
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
    index.mask_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&index.sa_transformed_mask);

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
    std::string merged = export_ms(a) + export_ms(b);
    if (a.k <= 32)  {
        return construct<uint64_t>(merged, a.k, a.klcp.size() > 0);
    } else {
        return construct<__uint128_t>(merged, a.k, a.klcp.size() > 0);
    }
}



void dump_index(const fms_index& index, const std::string &fn) {
    auto basename = fn + ".fmsi";
    sdsl::store_to_file(index.ac_gt, basename + ".ac_gt");
    sdsl::store_to_file(index.ac, basename + ".ac");
    sdsl::store_to_file(index.gt, basename + ".gt");
    sdsl::store_to_file(index.sa_transformed_mask, basename + ".mask");
    if (index.klcp.size() > 0) {
        sdsl::store_to_file(index.klcp, basename + ".klcp");
    }
    std::ofstream out(basename + ".misc");
    out << index.dollar_position << std::endl;
    for (auto c : index.counts) {
        out << c << std::endl;
    }
    out << index.k << std::endl;
    out.close();
}

fms_index load_index(const std::string &fn, bool use_klcp = true) {
    fms_index index;
    auto basename = fn + ".fmsi";
    sdsl::load_from_file(index.ac_gt, basename + ".ac_gt");
    index.ac_gt_rank = sdsl::rank_support_v5<1>(&index.ac_gt);
    sdsl::load_from_file(index.ac, basename + ".ac");
    index.ac_rank = sdsl::rank_support_v5<1>(&index.ac);
    sdsl::load_from_file(index.gt, basename + ".gt");
    index.gt_rank = sdsl::rank_support_v5<1>(&index.gt);
    sdsl::load_from_file(index.sa_transformed_mask, basename + ".mask");
    index.mask_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&index.sa_transformed_mask);
    if (std::filesystem::exists(basename + ".klcp") && use_klcp) {
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

// inline bool selects_loaded(const fms_index& index) {
//     return index.ac_gt_select != nullptr && index.ac_select != nullptr &&
//             index.gt_select != nullptr && index.mask_select != nullptr;
// }


// void init_selects(fms_index& index) {
//     if (index.ac_gt_select == nullptr) {
//         index.ac_gt_select = sdsl::select_support_hyb<1>(&index.ac_gt);
//     }
//     if (index.ac_select == nullptr) {
//         index.ac_select = sdsl::select_support_hyb<1>(&index.ac);
//     }
//     if (index.gt_select == nullptr) {
//         index.gt_select = sdsl::select_support_hyb<1>(&index.gt);
//     }
//     if (index.mask_select == nullptr) {
//         index.mask_select = sdsl::select_support_rrr<1>(&index.sa_transformed_mask);
//     }
// }
