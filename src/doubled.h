#pragma once

// The doubled text of FMSI 2 (index2 --fast): the masked superstring X, a separator, and the reverse
// complement Y of X. Every k-mer of the set then occurs in both orientations, so that a query is searched
// as it is, on one strand, and never as its reverse complement.
//
// The text is T = X 0 Y with n = |X|. The separator 0 sorts before the characters A, C, G, T (1 to 4) and
// after the end of T, so the suffix array of T is the generalized suffix array of {X, Y} in which a suffix
// of Y comes before an equal suffix of X, and no k-mer spans the two strings. The k-mer of Y at position j
// is the reverse complement of the one of X at n - j - k, and it is represented iff the latter is, so the
// mask of Y is the one of X reversed. The BWT has two sentinels, at the rows of the suffixes T[0..] and
// Y = T[n + 1..].
//
// Its MBWT (the BWT, the SA-transformed mask and the kLCP array) is computed from the rows of the suffix
// array by write_rows(), whichever way the suffix array was obtained.

#include <atomic>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "fms_index.h"
#include "parallel.h"

namespace doubled {

/// The codes of T: the separator 0 and A, C, G, T as 1 to 4, and MASK_BIT set at the represented k-mers
/// once the suffix arrays are constructed (set_mask). PADDING zero bytes follow T, for word-wise reading.
static const uint8_t MASK_BIT = 0x80, SYMBOL_BITS = 0x07;
static const size_t PADDING = 8;
/// No position (the row after the last one).
static const uint64_t NONE = ~(uint64_t)0;

struct text {
    /// Length of the superstring X.
    size_t n = 0;
    int k = 0;
    std::vector<uint8_t> codes;
    /// The mask of X (bit i set iff X[i] is upper case), kept until set_mask() moves it into the codes.
    sdsl::bit_vector mask;

    /// Length of T = X 0 Y.
    size_t length() const { return 2 * n + 1; }
};

/// The codes of T for the masked superstring; the mask is kept apart until set_mask().
text make_text(const std::string &ms, int k, unsigned threads = 1) {
    text t;
    t.n = ms.size();
    t.k = k;
    const size_t n = t.n;
    t.codes.assign(2 * n + 1 + PADDING, 0);
    t.mask = sdsl::bit_vector(n, 0);
    uint8_t *codes = t.codes.data();
    uint64_t *mask = t.mask.data();
    std::atomic<bool> invalid{false};
    // A chunk is a whole number of words of the mask.
    parallel_chunks(n, threads, parallel_chunk(64), [&](size_t begin, size_t end, unsigned) {
        for (size_t i = begin; i < end; ++i) {
            uint8_t code = 1 + nucleotideToInt[(uint8_t)ms[i]];
            if (code > 4) {
                invalid = true;
                code = 1;
            }
            codes[i] = code;
            // Y[j] is the complement of X[n - 1 - j] and at position n + 1 + j of T.
            codes[2 * n - i] = 5 - code;
            if (is_upper(ms[i])) mask[i >> 6] |= 1ULL << (i & 63);
        }
    });
    if (invalid) throw std::invalid_argument("the masked superstring contains a character other than ACGT");
    return t;
}

/// Mark the represented k-mers of X and Y in the codes; after the suffix arrays are constructed.
void set_mask(text &t, unsigned threads = 1) {
    const size_t n = t.n, k = t.k;
    uint8_t *codes = t.codes.data();
    const sdsl::bit_vector &mask = t.mask;
    parallel_chunks(n, threads, parallel_chunk(64), [&](size_t begin, size_t end, unsigned) {
        for (size_t i = begin; i < end; ++i) {
            if (!mask[i]) continue;
            codes[i] |= MASK_BIT;
            // The reverse complement of the k-mer of X at i is the k-mer of Y at n - i - k.
            if (i + k <= n) codes[2 * n + 1 - i - k] |= MASK_BIT;
        }
    });
    t.mask = sdsl::bit_vector();
}

/// Whether the suffix of T at p starts with k - 1 characters of X or of Y.
inline bool has_prefix(const text &t, uint64_t p) {
    const uint64_t m = t.k - 1;
    return p + m <= t.n || (p > t.n && p + m <= 2 * t.n + 1);
}

/// Whether the suffixes of T at p and q start with the same k - 1 symbols (the mask bits ignored).
inline bool same_prefix(const text &t, uint64_t p, uint64_t q) {
    const uint64_t symbols = 0x0707070707070707ULL;
    const uint8_t *a = t.codes.data() + p, *b = t.codes.data() + q;
    size_t m = t.k - 1, w = 0;
    for (; w + 8 <= m; w += 8) {
        uint64_t x, y;
        memcpy(&x, a + w, 8);
        memcpy(&y, b + w, 8);
        if ((x ^ y) & symbols) return false;
    }
    if (w == m) return true;
    uint64_t x, y;
    memcpy(&x, a + w, 8);
    memcpy(&y, b + w, 8);
    return ((x ^ y) & symbols & ((1ULL << (8 * (m - w))) - 1)) == 0;
}

/// The MBWT of T with all rows allocated and zero.
mbwt empty_mbwt(const text &t, bool use_klcp) {
    mbwt m;
    m.bwt.assign(t.length() + 1, 0);
    m.sa_transformed_mask = sdsl::bit_vector(t.length() + 1, 0);
    if (use_klcp) m.klcp = sdsl::bit_vector(t.length() + 1, 0);
    return m;
}

/// How far ahead the rows prefetch the text.
static const size_t ROW_PREFETCH = 32;

/// Fill the rows [row, row + count) of the MBWT, whose suffixes start at the given positions of T; [next]
/// is the position of the row after them (NONE if there is none), for the kLCP array. Rows of different
/// calls may share words of the bitvectors, so the words are written atomically.
template <typename position_t>
void write_rows(const text &t, mbwt &m, size_t row, const position_t *positions, size_t count, uint64_t next) {
    const uint8_t *codes = t.codes.data();
    const size_t n = t.n;
    byte *bwt = m.bwt.data();
    uint64_t *mask = m.sa_transformed_mask.data();
    uint64_t *klcp = m.klcp.size() > 0 ? m.klcp.data() : nullptr;
    uint64_t mask_word = 0, klcp_word = 0;
    size_t word = row >> 6;
    auto flush = [&]() {
        if (mask_word != 0) __atomic_fetch_or(mask + word, mask_word, __ATOMIC_RELAXED);
        if (klcp_word != 0) __atomic_fetch_or(klcp + word, klcp_word, __ATOMIC_RELAXED);
        mask_word = klcp_word = 0;
    };
    for (size_t x = 0; x < count; ++x) {
        if (x + ROW_PREFETCH < count) {
            uint64_t ahead = positions[x + ROW_PREFETCH];
            __builtin_prefetch(codes + ahead - (ahead > 0));
            if (klcp != nullptr) __builtin_prefetch(codes + ahead + t.k - 2);
        }
        size_t i = row + x;
        if ((i >> 6) != word) {
            flush();
            word = i >> 6;
        }
        uint64_t p = positions[x];
        uint8_t before = p == 0 ? 0 : codes[p - 1] & SYMBOL_BITS;
        if (before == 0) {
            // The suffix T[0..] or Y = T[n + 1..], preceded by the end of T or the separator.
            bwt[i] = 0;
            if (p == 0) {
                m.dollar_position = i;
            } else {
                m.dollar_position2 = i;
            }
        } else {
            bwt[i] = before - 1;
        }
        if (p < 2 * n + 1 && (codes[p] & MASK_BIT)) mask_word |= 1ULL << (i & 63);
        if (klcp != nullptr) {
            uint64_t q = x + 1 < count ? (uint64_t)positions[x + 1] : next;
            if (q != NONE && has_prefix(t, p) && has_prefix(t, q) && same_prefix(t, p, q)) {
                klcp_word |= 1ULL << (i & 63);
            }
        }
    }
    flush();
}

/// The MBWT of T from its suffix array, which includes the empty suffix (sa[0] = length of T).
template <typename sa_t>
void write_all_rows(const text &t, mbwt &m, const sa_t *sa, unsigned threads) {
    const size_t rows = t.length() + 1;
    parallel_chunks(rows, threads, parallel_chunk(64), [&](size_t begin, size_t end, unsigned) {
        write_rows(t, m, begin, sa + begin, end - begin, end < rows ? (uint64_t)sa[end] : NONE);
    });
}

/// Suffix array of the codes [0, length), including the empty suffix, by libsais.
template <typename sa_t>
std::vector<sa_t> suffix_array(const uint8_t *codes, size_t length, unsigned threads) {
    std::vector<sa_t> sa(length + 1);
    sa[0] = (sa_t)length;
    int64_t status;
#ifdef LIBSAIS_OPENMP
    if constexpr (sizeof(sa_t) == 4) {
        status = libsais_omp(codes, sa.data() + 1, (int32_t)length, 0, nullptr, (int32_t)threads);
    } else {
        status = libsais64_omp(codes, sa.data() + 1, (int64_t)length, 0, nullptr, (int64_t)threads);
    }
#else
    (void)threads;
    if constexpr (sizeof(sa_t) == 4) {
        status = libsais(codes, sa.data() + 1, (int32_t)length, 0, nullptr);
    } else {
        status = libsais64(codes, sa.data() + 1, (int64_t)length, 0, nullptr);
    }
#endif
    if (status != 0) throw std::runtime_error("libsais failed to construct the suffix array");
    return sa;
}

/// The MBWT of the doubled text of the masked superstring, which is consumed: the suffix array of T by
/// libsais, with 32-bit integers if T has less than 2^31 characters (5 bytes per character of T), else
/// with 64-bit ones (9 bytes per character).
mbwt compute_mbwt(std::string ms, int k, bool use_klcp, unsigned threads = 1) {
    text t = make_text(ms, k, threads);
    std::string().swap(ms);
    build_phases.lap("encode");
    mbwt m = empty_mbwt(t, use_klcp);
    if (t.length() < (size_t)INT32_MAX) {
        auto sa = suffix_array<int32_t>(t.codes.data(), t.length(), threads);
        build_phases.lap("suffix_array");
        set_mask(t, threads);
        write_all_rows(t, m, sa.data(), threads);
    } else {
        auto sa = suffix_array<int64_t>(t.codes.data(), t.length(), threads);
        build_phases.lap("suffix_array");
        set_mask(t, threads);
        write_all_rows(t, m, sa.data(), threads);
    }
    build_phases.lap(use_klcp ? "bwt_mask_klcp" : "bwt_mask");
    return m;
}

} // namespace doubled
