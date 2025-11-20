/// This is the main entry point for FMSI to be ebmedded within other tools.
/// Relying on other internals of FMSI is discouraged as they might be subject to change.
/// If there is anything missing, contact the authors.

#include "fms_index.h"

/**
 * Constructs FMSI.
 * 
 * @param ms a masked superstring in masked-case (ACGT implies 1 in mask, acgt 0 in mask).
 * @param k the length of one k-mer
 * @param construct_klcp whether the construction should also create a kLCP array.
 *                       Uses more memory (1 bit per char) and time, but enables streamed queries.
 * @param construct_access_support whether the construction should also create support for k-mer access.
 *                                 Uses more memory (<=0.4 bits per char) 
 * @param T the integer size used for construction of kLCP, should be smallest integers with bitlength larger than k.
 *          Relevant only if [construct_klcp] is set to `true`.
 */
template <typename T>
fms_index fmsi_construct(std::string &ms, int k, bool construct_klcp, bool construct_access_support) {
    fms_index index = construct<T>(ms, k, construct_klcp);
    if (construct_access_support) init_selects(index);
    return index;
}

/**
 * Constructs a support for kmer access queries. Needs to be executed before running access queries
 * if [construct_access_support] was not set to `true` during fmsi_construct.
 * Increases space by at most 0.8 bits / superstring char.
 */
void fmsi_construct_access_support(fms_index& index) {
    init_selects(index);
}

/**
 * Removes a support for kmer access queries to save space.
 * Decreases space by at most 0.8 bits / superstring char.
 */
void fmsi_deconstruct_access_support(fms_index& index) {
    destroy_selects(index);
}


/**
 * Performs a single membership query on the given k-mer.
 * Works in uni-directional model. In bi-directional model,
 * needs to be queried with the reverse complement if result is `false`.
 * 
 * @param index the FMSI index structure
 * @param kmer the queried k-mer. Ensure the length is exactly k
 * @param k the length of the k-mer.
 * @param maximized_ones whether the underlying masked superstring has mask maximizing ones
 *                       if this is `true`, the queries are faster
 *                       if the masked superstring does not maximize ones, then it has to be set to `false`, otherwise results may be incorrect
 *                       it is recommended, especially for benchmarking to first optimize the masked superstring and set this to `true`.
 * 
 * @return true if k-mer is present in the index.
 */
template <bool maximized_ones>
bool fmsi_membership_single_query(fms_index& index, char* kmer, int k) {
    return single_query_or<maximized_ones>(index, kmer, k) == 1;
}


/**
 * Performs a single dictionary query on the given k-mer.
 * Works in uni-directional model. In bi-directional model,
 * needs to be queried with the reverse complement if result is `-1`.
 * Important: if minimality of the hash is needed, set [minimal] to `true` and ensure the mask minimizes the number of ones.
 * 
 * @param index the FMSI index structure
 * @param kmer the queried k-mer. Ensure the length is exactly k
 * @param k the length of the k-mer.
 * @param minimal_hash whether the returned hashes should be minimal, i.e. from {0, |K|-1}
 *                     if minimality is not required, setting this to `false` increases speed
 * @param maximized_ones relevant only if [minimal_hash] is set to `false`
 *                       whether the underlying masked superstring has mask maximizing ones
 *                       if this is `true`, the queries are faster
 *                       if the masked superstring does not maximize ones, then it has to be set to `false`, otherwise results may be incorrect
 *                       it is recommended, especially for benchmarking to first optimize the masked superstring and set this to `true`.
 * 
 * @return if [minimal_hash] is `true`, for present k-mers returns an unique minimal hash from {0, |K|-1}.
 *         otherwise, for present k-mers returns an unique non-minimal hash from {0, |S|-1}.
 *         for absent k-mers, `-1` is always returned.
 */
template <bool minimal_hash, bool maximized_ones>
int64_t fmsi_lookup_single_query(fms_index& index, char* kmer, int k) {
    if constexpr (minimal_hash) {
        if constexpr (maximized_ones) {
            throw new std::invalid_argument("Minimal hashes are not possible on masked superstrings maximizing the number of ones."); 
        } else {
            return single_query_order(index, kmer, k);
        }
    } else {
        return single_query_order_nonminimal<maximized_ones>(index, kmer, k);
    }
}

/**
 * Performs a membership query for every k-mer in the pattern.
 * Works in bi-directional model.
 * Non-ACGT k-mers are treated as absent.
 * 
 * @param index the FMSI index structure. For the best performance, ensure kLCP array is constructed.
 * @param pattern the queried sequence
 * @param reverse_complementary_pattern should be the reverse complement of pattern
 * @param sequence_length the length of the queried pattern
 * @param k the length of the k-mer.
 * @param maximized_ones whether the underlying masked superstring has mask maximizing ones
 *                       if this is `true`, the queries are faster
 *                       if the masked superstring does not maximize ones, then it has to be set to `false`, otherwise results may be incorrect
 *                       it is recommended, especially for benchmarking to first optimize the masked superstring and set this to `true`.
 * 
 * @return a vector with 1 if the correcsponding k-mer is present in the index and 0 otherwise.
 */
template <bool maximized_ones>
std::vector<int64_t> fmsi_membership_streamed_query(fms_index& index, char* pattern, char* reverse_complementary_pattern, size_t sequence_length, int k) {
    return query_kmers_streaming_with_chunking<false, maximized_ones>(index, pattern, reverse_complementary_pattern, sequence_length, k, false);
}


/**
 * Performs a dictionary query for every k-mer in the pattern.
 * Works in bi-directional model.
 * Non-ACGT k-mers are treated as absent.
 * 
 * @param index the FMSI index structure. For the best performance, ensure kLCP array is constructed.
 * @param pattern the queried sequence
 * @param reverse_complementary_pattern should be the reverse complement of pattern
 * @param sequence_length the length of the queried pattern
 * @param k the length of the k-mer.
 * @param minimal_hash whether the returned hashes should be minimal, i.e. from {0, |K|-1}
 *                     if minimality is not required, setting this to `false` increases speed
 * @param maximized_ones relevant only if [minimal_hash] is set to `false`
 *                       whether the underlying masked superstring has mask maximizing ones
 *                       if this is `true`, the queries are faster
 *                       if the masked superstring does not maximize ones, then it has to be set to `false`, otherwise results may be incorrect
 *                       it is recommended, especially for benchmarking to first optimize the masked superstring and set this to `true`.
 * 
 * @return A vector with the following for each k-mer:
 *         if [minimal_hash] is `true`, for present k-mers returns an unique minimal hash from {0, |K|-1}.
 *         otherwise, for present k-mers returns an unique non-minimal hash from {0, |S|-1}.
 *         for absent k-mers, `-1` is always returned.
 */
template <bool minimal_hash, bool maximized_ones>
std::vector<int64_t> fmsi_lookup_streamed_query(fms_index& index, char* pattern, char* reverse_complementary_pattern, size_t sequence_length, int k) {
    return query_kmers_streaming_with_chunking<minimal_hash, maximized_ones>(index, pattern, reverse_complementary_pattern, sequence_length, k, true);
}

/**
 * Obtain a k-mer corresponding to a k-mer hash. The inverse of lookup.
 * Important: the minimal_hash parameter must be the same as for access.
 * 
 * @param index the FMSI index structure. Needs to have selects initialized.
 * @param hash the hash of k-mer. Ensure that this is non-negative.
 * @param k the length of k-mer.
 * @param minimal_hash true if the hash used was minimal (from 0 to |K|-1) or non-minimal (from 0 to |S|-1).
 * 
 * @return the k-mer corresponding to the hash.
 */
template <bool minimal_hash>
std::string fmsi_access(fms_index& index, int64_t hash, int k) {
    return kmer_access<minimal_hash>(index, hash, k);
}