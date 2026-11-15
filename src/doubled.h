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
// array by write_rows(), whichever way the suffix array was obtained: by libsais on T (compute_mbwt), or
// by merging the suffix arrays of X and of Y (compute_mbwt_lean), which need half the memory and 32-bit
// integers up to |X| < 2^31. Both give the same MBWT.

#include <atomic>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "fms_index.h"
#include "parallel.h"
#include "rank2.h"

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

/// Whether the suffix of T at p starts with m characters of X or of Y.
inline bool has_prefix(const text &t, uint64_t p, uint64_t m) {
    return p + m <= t.n || (p > t.n && p + m <= 2 * t.n + 1);
}

/// Whether the suffix of T at p starts with k - 1 characters of X or of Y.
inline bool has_prefix(const text &t, uint64_t p) { return has_prefix(t, p, t.k - 1); }

/// Whether the suffixes of T at p and q start with the same m symbols (the mask bits ignored).
inline bool same_prefix(const text &t, uint64_t p, uint64_t q, size_t m) {
    const uint64_t symbols = 0x0707070707070707ULL;
    const uint8_t *a = t.codes.data() + p, *b = t.codes.data() + q;
    size_t w = 0;
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

/// Whether the suffixes of T at p and q start with the same k - 1 symbols.
inline bool same_prefix(const text &t, uint64_t p, uint64_t q) { return same_prefix(t, p, q, t.k - 1); }

/// Whether the k-mer of T at p is its own reverse complement (the codes of complementary characters add up to 5).
inline bool is_palindrome(const text &t, uint64_t p) {
    const uint8_t *codes = t.codes.data() + p;
    for (int j = 0; j < t.k / 2; ++j) {
        if ((codes[j] & SYMBOL_BITS) + (codes[t.k - 1 - j] & SYMBOL_BITS) != 5) return false;
    }
    return t.k % 2 == 0;
}

/// The MBWT of T with all rows allocated and zero.
mbwt empty_mbwt(const text &t, bool use_klcp) {
    mbwt m;
    m.bwt.assign(t.length() + 1, 0);
    m.sa_transformed_mask = sdsl::bit_vector(t.length() + 1, 0);
    if (use_klcp) m.klcp = sdsl::bit_vector(t.length() + 1, 0);
    if (mbwt_with_same_kmer) m.same_kmer = sdsl::bit_vector(t.length() + 1, 0);
    if (mbwt_with_same_kmer && t.k % 2 == 0) m.palindrome = sdsl::bit_vector(t.length() + 1, 0);
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
    uint64_t *same = m.same_kmer.size() > 0 ? m.same_kmer.data() : nullptr;
    uint64_t *palindrome = m.palindrome.size() > 0 ? m.palindrome.data() : nullptr;
    uint64_t mask_word = 0, klcp_word = 0, same_word = 0, palindrome_word = 0;
    size_t word = row >> 6;
    auto flush = [&]() {
        if (mask_word != 0) __atomic_fetch_or(mask + word, mask_word, __ATOMIC_RELAXED);
        if (klcp_word != 0) __atomic_fetch_or(klcp + word, klcp_word, __ATOMIC_RELAXED);
        if (same_word != 0) __atomic_fetch_or(same + word, same_word, __ATOMIC_RELAXED);
        if (palindrome_word != 0) __atomic_fetch_or(palindrome + word, palindrome_word, __ATOMIC_RELAXED);
        mask_word = klcp_word = same_word = palindrome_word = 0;
    };
    for (size_t x = 0; x < count; ++x) {
        if (x + ROW_PREFETCH < count) {
            uint64_t ahead = positions[x + ROW_PREFETCH];
            __builtin_prefetch(codes + ahead - (ahead > 0));
            if (klcp != nullptr || same != nullptr) __builtin_prefetch(codes + ahead + t.k - 2);
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
        if (palindrome != nullptr && has_prefix(t, p, t.k) && is_palindrome(t, p)) palindrome_word |= 1ULL << (i & 63);
        if (klcp != nullptr || same != nullptr) {
            uint64_t q = x + 1 < count ? (uint64_t)positions[x + 1] : next;
            if (klcp != nullptr && q != NONE && has_prefix(t, p) && has_prefix(t, q) && same_prefix(t, p, q)) {
                klcp_word |= 1ULL << (i & 63);
            }
            if (same != nullptr && q != NONE && has_prefix(t, p, t.k) && has_prefix(t, q, t.k) && same_prefix(t, p, q, t.k)) {
                same_word |= 1ULL << (i & 63);
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
    if (t.length() < suffix_array_32bit_limit) {
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

/// LF-mapping over the BWT of X alone (its sentinel stored as A), for inserting the suffixes of Y.
struct lf_mapping {
    fmsi2::quad_sequence bwt;
    uint64_t counts[4] = {0, 0, 0, 0};
    uint64_t dollar = 0;

    /// The number of suffixes of X (with its empty one, which stands for the suffix of T at the separator)
    /// smaller than c S, given the number r of those smaller than S, for a suffix S of Y.
    inline uint64_t step(uint64_t r, uint8_t c) const {
        return counts[c] + bwt.rank(r, c) - (c == 0 && r > dollar);
    }
};

/// The rows of a suffix array of X or Y: sa[r] is the position, in X or Y, of the suffix of row r.
template <typename sa_t>
lf_mapping make_lf_mapping(const text &t, const sa_t *sa, unsigned threads) {
    const size_t rows = t.n + 1;
    const uint8_t *codes = t.codes.data();
    lf_mapping lf;
    lf.bwt = fmsi2::quad_sequence(rows, 64);
    std::vector<std::array<uint64_t, 4>> counted(threads, std::array<uint64_t, 4>{0, 0, 0, 0});
    std::vector<uint64_t> dollars(threads, NONE);
    // Four symbols make a byte, so the threads write to different bytes.
    parallel_chunks(rows, threads, parallel_chunk(4), [&](size_t begin, size_t end, unsigned thread) {
        for (size_t r = begin; r < end; ++r) {
            uint64_t p = sa[r];
            uint8_t c = p == 0 ? 0 : codes[p - 1] - 1;
            if (p == 0) dollars[thread] = r;
            lf.bwt.set(r, c);
            counted[thread][c]++;
        }
    });
    lf.bwt.finalize(threads);
    uint64_t totals[4] = {0, 0, 0, 0};
    for (unsigned thread = 0; thread < threads; ++thread) {
        for (int c = 0; c < 4; ++c) totals[c] += counted[thread][c];
        if (dollars[thread] != NONE) lf.dollar = dollars[thread];
    }
    lf.counts[0] = 1;
    lf.counts[1] = totals[0];
    lf.counts[2] = totals[0] + totals[1];
    lf.counts[3] = totals[0] + totals[1] + totals[2];
    return lf;
}

/// Whether the suffix of T at p (in X, or p = n for the one at the separator) is smaller than the
/// suffix of Y at j; the suffix of X ends with the separator at the latest, which Y does not contain.
inline bool x_suffix_smaller(const text &t, uint64_t p, uint64_t j) {
    const uint8_t *codes = t.codes.data();
    const uint64_t end = 2 * t.n + 1;
    for (uint64_t a = p, b = t.n + 1 + j;; ++a, ++b) {
        // The end of T sorts before the separator.
        if (b == end) return false;
        if (codes[a] != codes[b]) return codes[a] < codes[b];
    }
}

/// Write the array to a file, or read [count] items of it from [first] on.
template <typename T>
void write_array(const std::string &path, const T *items, size_t count) {
    FILE *file = fopen(path.c_str(), "wb");
    if (file == nullptr) throw std::runtime_error("couldn't write file " + path);
    size_t written = fwrite(items, sizeof(T), count, file);
    if (fclose(file) != 0 || written != count) throw std::runtime_error("couldn't write file " + path);
}

template <typename T>
bool read_array(int descriptor, size_t first, size_t count, T *items) {
    size_t bytes = count * sizeof(T), done = 0;
    char *out = reinterpret_cast<char *>(items);
    while (done < bytes) {
        ssize_t got = pread(descriptor, out + done, bytes - done, (off_t)(first * sizeof(T) + done));
        if (got <= 0) return false;
        done += (size_t)got;
    }
    return true;
}

/// Number of chains of the insertion of the suffixes of Y per thread, and how many of them a thread
/// advances at the same time.
static const size_t CHAINS_PER_THREAD = 64, INTERLEAVED_CHAINS = 16;

/// The MBWT of the doubled text of the masked superstring, which is consumed, by merging the suffix arrays
/// of X and Y (libsais with 32-bit integers, |X| < 2^31): T and one suffix array of n + 1 integers are in
/// memory at a time, and the suffix arrays are kept in temporary files ([temporary] followed by .sa_x and
/// .sa_y) until they are merged.
///
/// The suffixes of Y are inserted among those of X by backward search on the BWT of X: if r suffixes of X
/// are smaller than the suffix S of Y, then counts[c] + rank_c(r) are smaller than c S. Chains of these
/// steps run in parallel from anchors, suffixes of Y whose number is found by binary search on the suffix
/// array of X; each chain must arrive at the number of the next anchor, which checks the construction.
/// The numbers give, for each row of X, how many rows of Y come before it (the gaps), and the merged
/// order of the suffixes of T is then that of the generalized suffix array computed by compute_mbwt.
mbwt compute_mbwt_lean(std::string ms, int k, bool use_klcp, unsigned threads, const std::string &temporary) {
    if (ms.size() >= (size_t)INT32_MAX) {
        throw std::invalid_argument("the lean construction needs a superstring of less than 2^31 characters");
    }
    text t = make_text(ms, k, threads);
    std::string().swap(ms);
    build_phases.lap("encode");
    const size_t n = t.n;
    const uint8_t *codes = t.codes.data();
    const std::string path_x = temporary + ".sa_x", path_y = temporary + ".sa_y";

    // The suffix array of X, its LF-mapping and the numbers of the anchors.
    const size_t chains = std::max<size_t>(1, std::min<size_t>(n, (size_t)threads * CHAINS_PER_THREAD));
    // Anchor s is the suffix of Y at anchor_position(s); the last one is the empty suffix.
    auto anchor_position = [&](size_t s) { return (uint64_t)((unsigned __int128)s * n / chains); };
    std::vector<uint64_t> anchor_rank(chains + 1, 0);
    lf_mapping lf;
    {
        auto sa = suffix_array<int32_t>(codes, n, threads);
        build_phases.lap("suffix_array_x");
        lf = make_lf_mapping(t, sa.data(), threads);
        parallel_chunks(chains, threads, 1, [&](size_t begin, size_t end, unsigned) {
            for (size_t s = begin; s < end; ++s) {
                uint64_t j = anchor_position(s);
                anchor_rank[s] = std::partition_point(sa.begin(), sa.end(), [&](int32_t p) {
                    return x_suffix_smaller(t, (uint64_t)p, j);
                }) - sa.begin();
            }
        });
        write_array(path_x, sa.data(), sa.size());
        build_phases.lap("lf_anchors_x");
    }
    {
        auto sa = suffix_array<int32_t>(codes + n + 1, n, threads);
        build_phases.lap("suffix_array_y");
        write_array(path_y, sa.data(), sa.size());
        build_phases.lap("write_y");
    }

    // The gaps: gap[r] rows of Y come before the row r of X (r = n + 1: after the last one). They are
    // counted in 16 bits; the counts beyond go to the overflow of the thread.
    std::vector<uint16_t> gap(n + 2, 0);
    std::vector<std::vector<uint64_t>> overflow(threads);
    std::atomic<bool> inconsistent{false};
    // The empty suffix of Y is the smallest suffix of T.
    gap[0] = 1;
    auto count_gap = [&](uint64_t r, unsigned thread) {
        uint16_t value = __atomic_load_n(&gap[r], __ATOMIC_RELAXED);
        while (value != UINT16_MAX) {
            if (__atomic_compare_exchange_n(&gap[r], &value, (uint16_t)(value + 1), true, __ATOMIC_RELAXED,
                                            __ATOMIC_RELAXED)) {
                return;
            }
        }
        overflow[thread].push_back(r);
    };
    // Chain s inserts the suffixes of Y at the positions [anchor_position(s), anchor_position(s + 1)),
    // from the right, starting from the number of anchor s + 1.
    parallel_chunks(chains, threads, INTERLEAVED_CHAINS, [&](size_t begin, size_t end, unsigned thread) {
        struct chain {
            uint64_t r;
            int64_t j, stop;
            bool pending;
        };
        std::vector<chain> running;
        for (size_t s = begin; s < end; ++s) {
            running.push_back({anchor_rank[s + 1], (int64_t)anchor_position(s + 1) - 1, (int64_t)anchor_position(s),
                               false});
        }
        for (size_t active = running.size(); active > 0;) {
            active = 0;
            for (chain &c : running) {
                if (c.pending) {
                    // The gap and the block of the BWT were prefetched in the previous round.
                    count_gap(c.r, thread);
                    c.pending = false;
                }
                if (c.j < c.stop) continue;
                c.r = lf.step(c.r, codes[n + 1 + c.j] - 1);
                --c.j;
                c.pending = true;
                lf.bwt.prefetch(c.r);
                __builtin_prefetch(&gap[c.r], 1);
                ++active;
            }
        }
        for (size_t s = begin; s < end; ++s) {
            if (running[s - begin].r != anchor_rank[s]) inconsistent = true;
        }
    });
    if (inconsistent) throw std::runtime_error("the lean construction is inconsistent (chains and anchors differ)");
    std::unordered_map<uint64_t, uint64_t> extra;
    for (auto &list : overflow) {
        for (uint64_t r : list) extra[r]++;
        std::vector<uint64_t>().swap(list);
    }
    auto gap_of = [&](uint64_t r) -> uint64_t {
        return gap[r] == UINT16_MAX && extra.count(r) ? UINT16_MAX + extra.at(r) : gap[r];
    };
    lf = lf_mapping();
    build_phases.lap("insert_y");

    // The merge: the gaps are split into chunks, each knowing from their sums where its rows of Y and of T
    // start, and the rows of a chunk are written from the suffix arrays read from the files.
    set_mask(t, threads);
    mbwt m = empty_mbwt(t, use_klcp);
    const size_t gaps = n + 2, chunk = parallel_chunk_items, chunk_count = (gaps + chunk - 1) / chunk;
    std::vector<uint64_t> y_before(chunk_count + 1, 0);
    parallel_chunks(chunk_count, threads, 1, [&](size_t begin, size_t end, unsigned) {
        for (size_t c = begin; c < end; ++c) {
            uint64_t sum = 0;
            for (size_t r = c * chunk; r < std::min(gaps, (c + 1) * chunk); ++r) sum += gap_of(r);
            y_before[c + 1] = sum;
        }
    });
    for (size_t c = 0; c < chunk_count; ++c) y_before[c + 1] += y_before[c];
    if (y_before[chunk_count] != n + 1) throw std::runtime_error("the lean construction is inconsistent (gaps)");
    const size_t rows = t.length() + 1;
    int file_x = open(path_x.c_str(), O_RDONLY), file_y = open(path_y.c_str(), O_RDONLY);
    if (file_x < 0 || file_y < 0) throw std::runtime_error("couldn't read the temporary suffix arrays");
    std::atomic<bool> failed{false};
    parallel_chunks(chunk_count, threads, 1, [&](size_t begin, size_t end, unsigned) {
        std::vector<int32_t> xs, ys;
        std::vector<uint64_t> positions;
        for (size_t c = begin; c < end; ++c) {
            const size_t first = c * chunk, last = std::min(gaps, (c + 1) * chunk);
            const size_t x_first = first, x_last = std::min(last, n + 1);
            const size_t y_first = y_before[c], y_count = y_before[c + 1] - y_before[c];
            // One more suffix of each, for the row after the chunk.
            xs.resize(x_last - x_first + 1);
            ys.resize(y_count + 1);
            size_t x_read = std::min(x_last + 1, n + 1) - x_first, y_read = std::min(y_count + 1, n + 1 - y_first);
            if (!read_array(file_x, x_first, x_read, xs.data()) || !read_array(file_y, y_first, y_read, ys.data())) {
                failed = true;
                return;
            }
            positions.clear();
            size_t y = 0;
            for (size_t r = first; r < last; ++r) {
                for (uint64_t g = gap_of(r); g > 0; --g) positions.push_back(n + 1 + (uint64_t)ys[y++]);
                if (r <= n) positions.push_back((uint64_t)xs[r - x_first]);
            }
            const size_t row = first + y_first;
            uint64_t next = NONE;
            if (row + positions.size() < rows) {
                next = last < gaps && gap_of(last) > 0 ? n + 1 + (uint64_t)ys[y] : (uint64_t)xs[last - x_first];
            }
            write_rows(t, m, row, positions.data(), positions.size(), next);
        }
    });
    close(file_x);
    close(file_y);
    std::filesystem::remove(path_x);
    std::filesystem::remove(path_y);
    if (failed) throw std::runtime_error("couldn't read the temporary suffix arrays");
    build_phases.lap(use_klcp ? "merge_bwt_mask_klcp" : "merge_bwt_mask");
    return m;
}

} // namespace doubled
