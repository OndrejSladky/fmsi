#pragma once

// FMSI 2 prototype: the masked Burrows-Wheeler transform in a cache-friendly
// layout, searched in batches with prefetching.
//
// Differences from the original index (fms_index.h):
//  - the BWT is stored as 2-bit symbols with in-line rank counts (one cache
//    line per rank instead of a wavelet tree over three bitvectors),
//  - the SA-transformed mask is a plain bitvector with in-line rank counts, a
//    Tunstall-coded bitvector (smaller, one cache line per access as well),
//    or the original RRR vector,
//  - the SA intervals of all short patterns are tabulated when the index is loaded,
//  - k-mers are searched in lockstep, each step prefetching the cache lines of the next one.

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <sdsl/bit_vectors.hpp>

#include "fms_index.h"
#include "rank2.h"

namespace fmsi2 {

enum class mask_kind : uint32_t {
    plain = 0,
    rrr = 1,
    packed = 2,
};

struct index {
    uint32_t k = 0;
    /// Length of the BWT, i.e., the length of the superstring plus one.
    uint64_t size = 0;
    uint64_t dollar_position = 0;
    uint64_t counts[4] = {0, 0, 0, 0};
    quad_rank bwt;
    mask_kind mask_type = mask_kind::plain;
    bit_rank mask;
    tunstall_bits mask_packed;
    sdsl::rrr_vector<RRR_BLOCK_SIZE> mask_rrr;
    sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE> mask_rrr_rank;

    /// SA intervals of all patterns of length prefix_length; computed at load time, not stored.
    uint32_t prefix_length = 0;
    std::vector<uint64_t> prefix_table;

    index() = default;
    index(const index &) = delete;
    index &operator=(const index &) = delete;

    inline bool mask_get(size_t i) const {
        switch (mask_type) {
        case mask_kind::plain:
            return mask[i];
        case mask_kind::packed:
            return mask_packed[i];
        default:
            return mask_rrr[i];
        }
    }

    inline size_t mask_rank(size_t i) const {
        switch (mask_type) {
        case mask_kind::plain:
            return mask.rank(i);
        case mask_kind::packed:
            return mask_packed.rank(i);
        default:
            return mask_rrr_rank(i);
        }
    }

    /// Bit i, given the value returned by mask_prefetch(i).
    inline bool mask_get(size_t i, size_t hint) const {
        return mask_type == mask_kind::packed ? mask_packed.get(i, hint) : mask_get(i);
    }

    /// Prefetch what mask_get(i, hint) reads and return the hint.
    inline size_t mask_prefetch(size_t i) const {
        if (mask_type == mask_kind::plain) {
            mask.prefetch(i);
        } else if (mask_type == mask_kind::packed) {
            return mask_packed.prefetch(i);
        }
        return 0;
    }

    /// Whether mask_prefetch() reads memory that should be prefetched itself, by mask_prefetch_directory().
    inline bool mask_has_directory() const { return mask_type == mask_kind::packed; }

    inline void mask_prefetch_directory(size_t i) const {
        if (mask_type == mask_kind::packed) mask_packed.prefetch_directory(i);
    }

    size_t mask_bytes() const {
        switch (mask_type) {
        case mask_kind::plain:
            return mask.bytes();
        case mask_kind::packed:
            return mask_packed.bytes();
        default:
            return sdsl::size_in_bytes(mask_rrr);
        }
    }
};

/// Go from range (i,j) for pattern P to range for c+P.
inline void update_range(const index &ix, size_t &i, size_t &j, byte c) {
    size_t count = ix.counts[c];
    size_t rank_i, rank_j;
    ix.bwt.rank_pair(i, j, c, rank_i, rank_j);
    // The dollar is stored as A in the BWT.
    i = count + rank_i - (c == 0 && i > ix.dollar_position);
    j = count + rank_j - (c == 0 && j > ix.dollar_position);
}

/// Build the index from the MBWT, which is consumed.
void build(index &ix, mbwt &m, int k, mask_kind mask_type) {
    ix.k = k;
    ix.size = m.bwt.size();
    ix.dollar_position = m.dollar_position;

    ix.bwt = quad_rank(m.bwt.size());
    uint64_t totals[4] = {0, 0, 0, 0};
    for (size_t i = 0; i < m.bwt.size(); ++i) {
        ix.bwt.set(i, m.bwt[i]);
        totals[m.bwt[i]]++;
    }
    ix.bwt.finalize();
    m.bwt = std::vector<byte>();
    // The count of A includes the dollar, which stands for the smallest suffix.
    ix.counts[0] = 1;
    ix.counts[1] = totals[0];
    ix.counts[2] = totals[0] + totals[1];
    ix.counts[3] = totals[0] + totals[1] + totals[2];

    ix.mask_type = mask_type;
    if (mask_type == mask_kind::plain) {
        ix.mask = bit_rank(m.sa_transformed_mask.size());
        for (size_t i = 0; i < m.sa_transformed_mask.size(); ++i) {
            if (m.sa_transformed_mask[i]) ix.mask.set(i);
        }
        ix.mask.finalize();
    } else if (mask_type == mask_kind::packed) {
        ix.mask_packed.build(m.sa_transformed_mask.size(), [&](size_t i) { return (bool)m.sa_transformed_mask[i]; });
    } else {
        ix.mask_rrr = sdsl::rrr_vector<RRR_BLOCK_SIZE>(m.sa_transformed_mask);
        ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
    }
    m.sa_transformed_mask.resize(0);
}

static const char MAGIC[8] = {'F', 'M', 'S', 'I', '2', 'I', 'D', 'X'};
static const uint32_t FORMAT_VERSION = 1;

void save(const index &ix, const std::string &path) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::invalid_argument("couldn't open file " + path);
    out.write(MAGIC, sizeof(MAGIC));
    write_pod(out, FORMAT_VERSION);
    write_pod(out, ix.k);
    write_pod(out, ix.size);
    write_pod(out, ix.dollar_position);
    for (int c = 0; c < 4; ++c) write_pod(out, ix.counts[c]);
    write_pod(out, (uint32_t)ix.mask_type);
    ix.bwt.write(out);
    if (ix.mask_type == mask_kind::plain) {
        ix.mask.write(out);
    } else if (ix.mask_type == mask_kind::packed) {
        ix.mask_packed.write(out);
    } else {
        ix.mask_rrr.serialize(out);
    }
}

void load(index &ix, const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::invalid_argument("couldn't open file " + path);
    char magic[sizeof(MAGIC)];
    in.read(magic, sizeof(MAGIC));
    uint32_t version = 0, mask_type = 0;
    read_pod(in, version);
    if (!in || memcmp(magic, MAGIC, sizeof(MAGIC)) != 0 || version != FORMAT_VERSION) {
        throw std::invalid_argument("file " + path + " is not an FMSI 2 index of a supported version");
    }
    read_pod(in, ix.k);
    read_pod(in, ix.size);
    read_pod(in, ix.dollar_position);
    for (int c = 0; c < 4; ++c) read_pod(in, ix.counts[c]);
    read_pod(in, mask_type);
    ix.mask_type = (mask_kind)mask_type;
    ix.bwt.read(in);
    if (ix.mask_type == mask_kind::plain) {
        ix.mask.read(in);
    } else if (ix.mask_type == mask_kind::packed) {
        ix.mask_packed.read(in);
    } else {
        ix.mask_rrr.load(in);
        ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
    }
    if (!in) throw std::invalid_argument("file " + path + " is truncated");
}

/// Tabulate the SA intervals of all patterns of the given length. An entry is indexed by the
/// characters in the order of the backward search, the first one being the most significant.
void build_prefix_table(index &ix, unsigned length) {
    length = std::min(length, ix.k);
    std::vector<uint64_t> current = {0, ix.size};
    for (unsigned depth = 0; depth < length; ++depth) {
        std::vector<uint64_t> next(current.size() * 4);
        for (size_t key = 0; key < current.size() / 2; ++key) {
            for (byte c = 0; c < 4; ++c) {
                size_t i = current[2 * key], j = current[2 * key + 1];
                if (i != j) update_range(ix, i, j, c);
                next[2 * (4 * key + c)] = i;
                next[2 * (4 * key + c) + 1] = j;
            }
        }
        current = std::move(next);
    }
    ix.prefix_length = length;
    ix.prefix_table = std::move(current);
}

struct lookup_stats {
    uint64_t kmers = 0;
    uint64_t found = 0;
    /// Number of steps of the backward search (LF-mappings of an interval).
    uint64_t steps = 0;
    /// Number of searches started from scratch, and those of them caused by an unset mask bit of a long match.
    uint64_t searches = 0;
    uint64_t rechecks = 0;
    double seconds = 0;

    void add(const lookup_stats &other) {
        kmers += other.kmers;
        found += other.found;
        steps += other.steps;
        searches += other.searches;
        rechecks += other.rechecks;
    }
};

/// Look up the k-mers starting at the given offsets of the sequence, [batch_size] of them at a time.
/// For each k-mer, out[offset] is set to its presence (0/1), or to its hash (-1 if absent) if [hashes].
/// A k-mer is represented if the k-mer or its reverse complement has an occurrence set in the mask.
template <bool hashes, typename out_t>
void lookup_batch(const index &ix, const char *sequence, const uint32_t *offsets, size_t count, out_t *out,
                  unsigned batch_size, lookup_stats &stats) {
    const unsigned k = ix.k;
    const unsigned p = ix.prefix_length;
    const uint64_t *table = ix.prefix_table.data();

    enum : uint8_t { IDLE, START, SEARCH, LOCATE_MASK, MASK };
    struct slot {
        size_t i, j;
        /// The value returned by mask_prefetch(i).
        size_t hint;
        const char *kmer;
        uint64_t key;
        uint32_t offset;
        /// Number of characters of the k-mer already searched.
        uint32_t t;
        uint8_t state;
        /// Whether the reverse complement is being searched, and whether it is the second strand tried.
        uint8_t reverse, second;
    };
    std::vector<slot> slots(std::max(1u, batch_size));
    size_t next = 0, active = slots.size();
    // Saturating counter predicting the strand of the next k-mer from the previous hits.
    int score = 0;
    uint64_t found = 0, steps = 0, searches = 0;

    auto begin_strand = [&](slot &s) {
        ++searches;
        if (p == 0) {
            s.i = 0;
            s.j = ix.size;
            s.t = 0;
            s.state = SEARCH;
            return;
        }
        uint64_t key = 0;
        if (s.reverse) {
            for (unsigned t = 0; t < p; ++t) key = (key << 2) | (3 - nucleotideToInt[(uint8_t)s.kmer[t]]);
        } else {
            for (unsigned t = 0; t < p; ++t) key = (key << 2) | nucleotideToInt[(uint8_t)s.kmer[k - 1 - t]];
        }
        s.key = key;
        s.state = START;
        __builtin_prefetch(table + 2 * key);
    };
    auto load = [&](slot &s) {
        if (next == count) {
            s.state = IDLE;
            --active;
            return;
        }
        s.offset = offsets[next++];
        s.kmer = sequence + s.offset;
        s.reverse = score < 0;
        s.second = 0;
        begin_strand(s);
    };
    auto strand_failed = [&](slot &s) {
        if (s.second) {
            out[s.offset] = hashes ? -1 : 0;
            load(s);
        } else {
            s.second = 1;
            s.reverse ^= 1;
            begin_strand(s);
        }
    };
    const bool mask_has_directory = ix.mask_has_directory();
    auto advance = [&](slot &s) {
        if (s.i == s.j) {
            strand_failed(s);
        } else if (s.t == k && mask_has_directory) {
            // The directory of the mask is hardly in the cache after a whole search.
            s.state = LOCATE_MASK;
            ix.mask_prefetch_directory(s.i);
        } else if (s.t == k) {
            s.state = MASK;
            s.hint = ix.mask_prefetch(s.i);
        } else {
            s.state = SEARCH;
            ix.bwt.prefetch(s.i);
            ix.bwt.prefetch(s.j);
        }
    };

    for (slot &s : slots) load(s);
    while (active > 0) {
        for (slot &s : slots) {
            switch (s.state) {
            case START:
                s.i = table[2 * s.key];
                s.j = table[2 * s.key + 1];
                s.t = p;
                advance(s);
                break;
            case SEARCH: {
                byte c = s.reverse ? 3 - nucleotideToInt[(uint8_t)s.kmer[s.t]]
                                   : nucleotideToInt[(uint8_t)s.kmer[k - 1 - s.t]];
                update_range(ix, s.i, s.j, c);
                ++s.t;
                ++steps;
                advance(s);
                break;
            }
            case LOCATE_MASK:
                s.state = MASK;
                s.hint = ix.mask_prefetch(s.i);
                break;
            case MASK: {
                bool represented = ix.mask_get(s.i, s.hint);
                for (size_t q = s.i + 1; q < s.j && !represented; ++q) represented = ix.mask_get(q);
                if (represented) {
                    out[s.offset] = hashes ? (out_t)ix.mask_rank(s.i) : (out_t)1;
                    ++found;
                    score = std::max(-7, std::min(7, score + (s.reverse ? -1 : 1)));
                    load(s);
                } else {
                    strand_failed(s);
                }
                break;
            }
            default:
                break;
            }
        }
    }

    stats.kmers += count;
    stats.found += found;
    stats.steps += steps;
    stats.searches += searches;
}

/// A stretch of a query sequence without invalid characters, at least k characters long.
struct segment {
    uint32_t offset, length;
};

/// Look up all k-mers of the segments by streaming, [batch_size] segments at a time; the results
/// are as in lookup_batch (out must be initialized to 0, or to -1 if [hashes]).
///
/// A segment is searched on one strand and the k-mers not found there on the other one. On a strand,
/// the k-mers (windows) are resolved from the right to the left: the match of a window is extended by one
/// character to the match of the next one, so the match can get longer than k. A longer match has a
/// subset of the occurrences of the k-mer, hence a set mask bit proves the k-mer represented, while an
/// unset one requires searching the k-mer itself. This needs no kLCP array.
///
/// When a window cannot be found, the search stopped on a substring that does not occur in the
/// superstring, so that all windows containing it are absent on this strand. To find such substrings
/// cheaply, the unresolved windows containing the position where the search stopped are searched
/// starting from the lowest one, whose search stops the earliest.
template <bool hashes, typename out_t>
void lookup_streaming(const index &ix, const char *sequence, const segment *segments, size_t count, out_t *out,
                      unsigned batch_size, lookup_stats &stats) {
    const int k = ix.k;
    const int p = ix.prefix_length;
    const uint64_t *table = ix.prefix_table.data();
    // Values of out while a segment is processed; a non-negative hash or 1 means found.
    const out_t UNRESOLVED = hashes ? -1 : 0;
    const out_t ABSENT[2] = {(out_t)(hashes ? -2 : 2), (out_t)(hashes ? -3 : 3)};
    const int NO_BLOCK = -1;

    enum : uint8_t { IDLE, START, SEARCH, CHECK };
    struct slot {
        size_t i, j;
        /// The value returned by mask_prefetch(i).
        size_t hint;
        const char *sequence;
        /// Results of the windows of the segment.
        out_t *out;
        uint64_t key;
        int32_t length, windows;
        /// The window searched or matched, in the coordinates of the strand.
        int32_t x;
        /// Number of characters of the window already searched.
        int32_t t;
        /// The highest unresolved window if the lowest one of a block is searched, NO_BLOCK otherwise.
        int32_t block_top;
        /// Windows above this one are searched from the top only.
        int32_t top_only_above;
        int32_t hits[2];
        uint8_t state;
        /// Whether the reverse complement is searched, and whether it is the second strand,
        /// and whether the matched string has exactly k characters.
        uint8_t reverse, second, exact;
    };
    std::vector<slot> slots(std::max(1u, batch_size));
    size_t next = 0, active = slots.size();
    int score = 0;
    uint64_t found = 0, steps = 0, kmers = 0, searches = 0, rechecks = 0;

    // Character of the strand at position v, and the result of its window y.
    auto symbol = [&](const slot &s, int v) -> byte {
        return s.reverse ? 3 - nucleotideToInt[(uint8_t)s.sequence[s.length - 1 - v]]
                         : nucleotideToInt[(uint8_t)s.sequence[v]];
    };
    auto result = [&](const slot &s, int y) -> out_t & { return s.out[s.reverse ? s.windows - 1 - y : y]; };
    auto unresolved = [&](const slot &s, int y) {
        out_t value = result(s, y);
        return value == UNRESOLVED || (s.second && value == ABSENT[0]);
    };

    std::function<void(slot &)> load;

    auto begin_search = [&](slot &s, int y) {
        ++searches;
        s.x = y;
        s.t = 0;
        if (p == 0) {
            s.i = 0;
            s.j = ix.size;
            s.state = SEARCH;
            return;
        }
        uint64_t key = 0;
        for (int t = 0; t < p; ++t) key = (key << 2) | symbol(s, y + k - 1 - t);
        s.key = key;
        s.state = START;
        __builtin_prefetch(table + 2 * key);
    };
    std::function<void(slot &)> strand_done;
    // Search the highest unresolved window not above [from].
    auto next_top = [&](slot &s, int from) {
        while (from >= 0 && !unresolved(s, from)) --from;
        s.block_top = NO_BLOCK;
        if (from < 0) {
            strand_done(s);
        } else {
            begin_search(s, from);
        }
    };
    strand_done = [&](slot &s) {
        if (!s.second) {
            s.second = 1;
            s.reverse ^= 1;
            s.top_only_above = INT32_MAX;
            next_top(s, s.windows - 1);
        } else {
            for (int y = 0; y < s.windows; ++y) {
                if (s.out[y] == ABSENT[0] || s.out[y] == ABSENT[1]) s.out[y] = UNRESOLVED;
            }
            if (s.hits[0] != s.hits[1]) score = std::max(-7, std::min(7, score + (s.hits[0] > s.hits[1] ? 1 : -1)));
            load(s);
        }
    };
    load = [&](slot &s) {
        if (next == count) {
            s.state = IDLE;
            --active;
            return;
        }
        const segment &g = segments[next++];
        s.sequence = sequence + g.offset;
        s.out = out + g.offset;
        s.length = g.length;
        s.windows = g.length - k + 1;
        kmers += s.windows;
        s.reverse = score < 0;
        s.second = 0;
        s.hits[0] = s.hits[1] = 0;
        s.top_only_above = INT32_MAX;
        next_top(s, s.windows - 1);
    };
    // The search of window x stopped: the substring from position a to the end of the window does not occur.
    auto search_failed = [&](slot &s) {
        int a = s.x + k - s.t;
        if (s.block_top == NO_BLOCK) {
            result(s, s.x) = ABSENT[s.second];
            // The unresolved windows containing position a.
            int lowest = std::max(0, a - k + 1);
            int top = s.x - 1;
            while (top >= lowest && !unresolved(s, top)) --top;
            if (top < lowest) {
                next_top(s, top);
                return;
            }
            int bottom = lowest;
            while (!unresolved(s, bottom)) ++bottom;
            if (bottom == top || top > s.top_only_above) {
                begin_search(s, top);
            } else {
                s.block_top = top;
                begin_search(s, bottom);
            }
        } else {
            // All windows from x to a contain the substring.
            for (int y = s.x; y <= std::min(a, s.block_top); ++y) {
                if (unresolved(s, y)) result(s, y) = ABSENT[s.second];
            }
            int bottom = a + 1;
            while (bottom <= s.block_top && !unresolved(s, bottom)) ++bottom;
            if (bottom > s.block_top) {
                next_top(s, s.x - 1);
            } else if (bottom == s.block_top) {
                s.block_top = NO_BLOCK;
                begin_search(s, bottom);
            } else {
                begin_search(s, bottom);
            }
        }
    };
    auto to_check = [&](slot &s) {
        s.state = CHECK;
        s.hint = ix.mask_prefetch(s.i);
        ix.bwt.prefetch(s.i);
        ix.bwt.prefetch(s.j);
    };
    auto after_search_step = [&](slot &s) {
        if (s.i == s.j) {
            search_failed(s);
        } else if (s.t == k) {
            s.exact = 1;
            to_check(s);
        } else {
            s.state = SEARCH;
            ix.bwt.prefetch(s.i);
            ix.bwt.prefetch(s.j);
        }
    };
    // The interval is the one of a string starting with window x: resolve the window and go to the next one.
    auto check = [&](slot &s) {
        size_t q = s.i;
        if (!ix.mask_get(q, s.hint)) {
            ++q;
            while (q < s.j && !ix.mask_get(q)) ++q;
        }
        if (q < s.j) {
            result(s, s.x) = hashes ? (out_t)ix.mask_rank(q) : (out_t)1;
            ++found;
            ++s.hits[s.reverse];
        } else if (!s.exact) {
            // Other occurrences of the k-mer may be set in the mask.
            ++rechecks;
            s.block_top = NO_BLOCK;
            begin_search(s, s.x);
            return;
        } else {
            result(s, s.x) = ABSENT[s.second];
        }
        if (s.block_top != NO_BLOCK) {
            // The lowest window of a block occurs, so nothing is known about the windows above it.
            int top = s.block_top;
            s.top_only_above = s.x;
            s.block_top = NO_BLOCK;
            begin_search(s, top);
            return;
        }
        if (s.x == 0 || !unresolved(s, s.x - 1)) {
            next_top(s, s.x - 1);
            return;
        }
        --s.x;
        update_range(ix, s.i, s.j, symbol(s, s.x));
        ++steps;
        if (s.i == s.j) {
            begin_search(s, s.x);
        } else {
            s.exact = 0;
            to_check(s);
        }
    };

    for (slot &s : slots) load(s);
    while (active > 0) {
        for (slot &s : slots) {
            switch (s.state) {
            case START:
                s.i = table[2 * s.key];
                s.j = table[2 * s.key + 1];
                s.t = p;
                after_search_step(s);
                break;
            case SEARCH:
                update_range(ix, s.i, s.j, symbol(s, s.x + k - 1 - s.t));
                ++s.t;
                ++steps;
                after_search_step(s);
                break;
            case CHECK:
                check(s);
                break;
            default:
                break;
            }
        }
    }

    stats.kmers += kmers;
    stats.found += found;
    stats.steps += steps;
    stats.searches += searches;
    stats.rechecks += rechecks;
}

} // namespace fmsi2
