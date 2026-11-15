#pragma once

// FMSI 2 prototype: the masked Burrows-Wheeler transform in a cache-friendly
// layout, searched in batches with prefetching.
//
// Differences from the original index (fms_index.h):
//  - the BWT is stored as 2-bit symbols with in-line rank counts (one block of
//    64 to 256 bytes per rank instead of a wavelet tree over three bitvectors),
//  - the SA-transformed mask is a plain bitvector with in-line rank counts, a
//    Tunstall-coded bitvector (smaller, one cache line per access as well),
//    or the original RRR vector,
//  - the SA intervals of all short patterns are tabulated when the index is loaded,
//  - k-mers are searched in lockstep, each step prefetching the cache lines of the next one.

#include <array>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <sdsl/bit_vectors.hpp>

#include "fms_index.h"
#include "index_file.h"
#include "rank2.h"
#include "text_mask.h"
#include "version.h"

namespace fmsi2 {

enum class mask_kind : uint32_t {
    plain = 0,
    rrr = 1,
    packed = 2,
};

/// A position that is not a row.
static constexpr uint64_t NO_ROW = ~(uint64_t)0;

struct index {
    uint32_t k = 0;
    /// Number of strands of the superstring indexed: 1, or 2 if its reverse complement is indexed too
    /// (index2 --fast), in which case a k-mer is searched as it is and not as its reverse complement.
    uint32_t strands = 1;
    /// Length of the BWT, i.e., the number of suffixes including the empty one.
    uint64_t size = 0;
    /// Rows whose BWT symbol is a sentinel (stored as A): the one before the text and, if the text consists
    /// of two strings (the superstring and its reverse complement), the one between them; else NO_ROW.
    uint64_t dollar_position = 0;
    uint64_t dollar_position2 = NO_ROW;
    /// counts[c] is the number of suffixes starting with a symbol smaller than c, including those starting
    /// with a sentinel (the empty suffix and, with two strings, the one starting with the separator).
    uint64_t counts[4] = {0, 0, 0, 0};
    quad_sequence bwt;
    mask_kind mask_type = mask_kind::plain;
    bit_rank mask;
    tunstall_bits mask_packed;
    sdsl::rrr_vector<RRR_BLOCK_SIZE> mask_rrr;
    sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE> mask_rrr_rank;

    /// Optional kLCP array for streaming, stored in a file of its own: bit i is set iff the suffixes
    /// of the rows i and i + 1 share their first k - 1 characters.
    sdsl::bit_vector klcp;

    /// SA intervals of all patterns of length prefix_length, two integers per pattern (32-bit if the BWT has
    /// less than 2^32 rows); computed at load time, not stored.
    uint32_t prefix_length = 0;
    std::vector<uint32_t> prefix_table32;
    std::vector<uint64_t> prefix_table64;

    /// The interval of the pattern number key of the table.
    inline void prefix_interval(uint64_t key, size_t &i, size_t &j) const {
        if (!prefix_table32.empty()) {
            i = prefix_table32[2 * key];
            j = prefix_table32[2 * key + 1];
        } else {
            i = prefix_table64[2 * key];
            j = prefix_table64[2 * key + 1];
        }
    }

    inline void prefix_prefetch(uint64_t key) const {
        if (!prefix_table32.empty()) {
            __builtin_prefetch(prefix_table32.data() + 2 * key);
        } else {
            __builtin_prefetch(prefix_table64.data() + 2 * key);
        }
    }

    size_t prefix_table_bytes() const { return prefix_table32.size() * 4 + prefix_table64.size() * 8; }

    /// For access (hash -> k-mer), computed when the index is built or loaded, not stored: the number
    /// of symbols A of the BWT before the dollar, which is stored as A, and the number of set bits of
    /// the mask, i.e., of hashes.
    uint64_t dollar_rank = 0;
    uint64_t mask_ones = 0;
    /// The same for the second sentinel (NO_ROW if there is none); the sentinels ordered by their rows.
    uint64_t dollar_rank2 = NO_ROW;

    /// The preset the index was built with ("small" or "fast", "" if unknown) and the identifier of its content
    /// (0 if unknown), by which its kLCP file is recognized.
    std::string preset;
    uint64_t uid = 0;
    /// The type of the mask as the index file states it: "min-one" (every k-mer of the set has one set
    /// occurrence, so that lookup is a minimal perfect hash), "max-one" (all occurrences of the k-mers of the
    /// set are set), "min-one,max-one", "other", or "" if unknown; and how it is known ("computed", "declared").
    std::string mask_type_name, mask_type_source;
    /// The FASTA header of the masked superstring that was indexed, if the index file has it.
    std::string superstring_header;

    /// Bytes of the index in memory: the BWT, the mask, the kLCP array, the table of patterns and the hints.
    size_t bytes() const {
        return bwt.bytes() + mask_bytes() + (klcp.size() + 7) / 8 + prefix_table_bytes() + bwt.select_hint_bytes();
    }

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

    /// Position of the set bit number h (from 0) of the mask, for h < mask_ones.
    inline size_t mask_select(size_t h) const {
        if (mask_type == mask_kind::plain) return mask.select(h);
        // The other encodings have no select: the position is searched by rank from an interpolated guess.
        size_t guess = (size_t)((unsigned __int128)h * size / std::max<uint64_t>(1, mask_ones));
        return first_above(0, size, guess, [&](size_t i) { return mask_rank(i + 1) > h; });
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

    /// Whether a bit of the mask in [i, j) is set, for i < j, given the value returned by mask_prefetch(i).
    inline bool mask_any(size_t i, size_t j, size_t hint) const {
        if (mask_get(i, hint)) return true;
        // Most intervals are short; a long one (a k-mer with many occurrences) costs two ranks, not a scan.
        size_t stop = std::min(j, i + 4);
        for (size_t q = i + 1; q < stop; ++q) {
            if (mask_get(q)) return true;
        }
        return stop < j && mask_rank(stop) != mask_rank(j);
    }

    inline bool has_klcp() const { return klcp.size() > 0; }

    /// Whether the kLCP bits of rows i - 1 and i are unset, i.e., the first k - 1 characters of row i are those
    /// of no other row, so that klcp_widen leaves the interval [i, i + 1) as it is (i > 0).
    FMSI2_INLINE bool klcp_single(size_t i) const {
        const uint64_t *words = klcp.data();
        const size_t q = i - 1;
        return (((words[q >> 6] >> (q & 63)) | (words[i >> 6] >> (i & 63) << 1)) & 3) == 0;
    }

    FMSI2_INLINE void klcp_prefetch(size_t i, size_t j) const {
        __builtin_prefetch(klcp.data() + ((i - 1) >> 6));
        __builtin_prefetch(klcp.data() + ((j - 1) >> 6));
    }

    /// Go from the interval of a string of at least k - 1 characters to the one of its first k - 1 characters.
    FMSI2_INLINE void klcp_widen(size_t &i, size_t &j) const {
        const uint64_t *words = klcp.data();
        // The first unset bit from j - 1 on; the last bit of the array is unset.
        size_t w = (j - 1) >> 6;
        uint64_t unset = ~words[w] & (~0ULL << ((j - 1) & 63));
        while (unset == 0) unset = ~words[++w];
        j = (w << 6) + __builtin_ctzll(unset) + 1;
        // The last unset bit up to i - 1; bit 0 is unset.
        w = (i - 1) >> 6;
        unset = ~words[w] & (~0ULL >> (63 - ((i - 1) & 63)));
        while (unset == 0) unset = ~words[--w];
        i = (w << 6) + 64 - __builtin_clzll(unset);
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

/// Number of sentinels among the first i symbols of the BWT (they are stored as A).
inline size_t dollars_before(const index &ix, size_t i) {
    return (size_t)(i > ix.dollar_position) + (size_t)(i > ix.dollar_position2);
}

/// Go from range (i,j) for pattern P to range for c+P.
inline void update_range(const index &ix, size_t &i, size_t &j, byte c) {
    size_t count = ix.counts[c];
    size_t rank_i, rank_j;
    ix.bwt.rank_pair(i, j, c, rank_i, rank_j);
    // The sentinels are stored as A in the BWT.
    size_t dollars_i = c == 0 ? dollars_before(ix, i) : 0, dollars_j = c == 0 ? dollars_before(ix, j) : 0;
    i = count + rank_i - dollars_i;
    j = count + rank_j - dollars_j;
}

/// Set what access needs beyond the stored index and costs nothing.
inline void prepare_access(index &ix) {
    uint64_t first = std::min(ix.dollar_position, ix.dollar_position2);
    uint64_t second = std::max(ix.dollar_position, ix.dollar_position2);
    ix.dollar_rank = ix.bwt.rank(first, 0);
    ix.dollar_rank2 = second == NO_ROW ? NO_ROW : ix.bwt.rank(second, 0);
    ix.mask_ones = ix.mask_rank(ix.size);
}

/// Build the select hints of the BWT, which make access several times faster. They take one pass over
/// the blocks of the BWT and 1/16 bit per character, and are not stored. Without them, access works
/// too, guessing from the superblocks.
inline void build_access_support(index &ix) { ix.bwt.build_select_hints(); }

/// The type of a mask. FMSI answers membership queries the fastest on a max-one mask (all occurrences of the
/// k-mers of the set are set, so the first bit of an interval decides) and dictionary queries on a min-one mask
/// (one set occurrence per k-mer, so the rank of the interval is a minimal perfect hash). FMSI 2 tests any mask
/// in constant time, so the type matters for the hashes only: with a mask other than min-one they stay unique
/// per set occurrence, but a k-mer set several times has several, and their range exceeds the number of k-mers.
///
/// The type is computed from the rows that start with the same k-mer (mbwt::same_kmer). For an index of both
/// strands these hold the occurrences of a k-mer and of its reverse complement, so the type is exact. For an
/// index of one strand a k-mer and its reverse complement are not related to each other: a mask found min-one
/// may still set a k-mer once in each orientation (masks of KmerCamel in its bidirectional mode never do). The
/// type that the producer of the masked superstring states in its FASTA header (mask=min-one) or the user gives
/// is therefore preferred, and what is computed checks it.
struct mask_statistics {
    /// "min-one", "max-one", "min-one,max-one" (each k-mer of the set occurs once), "other" or "unknown".
    std::string type = "unknown";
    /// "computed" (exact: both strands), "computed-one-strand", "header" (of the masked superstring),
    /// "declared" (by the user), or empty.
    std::string source;
    /// Whether the type was computed, and what was found: no k-mer with several set occurrences; no k-mer of
    /// the set with an unset occurrence.
    bool computed = false, single = false, uniform = false;
    /// Number of different k-mers of the text with a set occurrence (computed only).
    uint64_t kmers = 0;
    /// Set if a stated type contradicts what was computed.
    std::string warning;
};

/// The type of the mask computed from the MBWT (needs mbwt::same_kmer; "unknown" without it).
inline mask_statistics classify_mask(const mbwt &m) {
    mask_statistics result;
    if (m.bwt.empty() || m.same_kmer.size() != m.bwt.size()) return result;
    bool single = true, uniform = true;
    uint64_t ones = 0, length = 0;
    const size_t rows = m.bwt.size();
    const bool palindromes = m.palindrome.size() == rows;
    for (size_t i = 0; i < rows; ++i) {
        ones += m.sa_transformed_mask[i];
        ++length;
        if (!m.same_kmer[i]) {
            if (ones > 0) {
                ++result.kmers;
                // A k-mer that is its own reverse complement occurs twice among the rows of both strands.
                const uint64_t once = palindromes && m.palindrome[i] ? 2 : 1;
                single = single && ones == once;
                uniform = uniform && ones == length;
            }
            ones = length = 0;
        }
    }
    result.computed = true;
    result.single = single;
    result.uniform = uniform;
    result.type = single && uniform ? "min-one,max-one" : single ? "min-one" : uniform ? "max-one" : "other";
    result.source = m.dollar_position2 == mbwt::NO_ROW ? "computed-one-strand" : "computed";
    return result;
}

/// The type of the mask given what was computed and what was stated: by the user ([declared]: "min-one",
/// "max-one", or "auto" for nothing) or in the header of the masked superstring ([header_type], may be empty).
/// A stated type is kept unless the computed one contradicts it.
inline mask_statistics decide_mask_type(mask_statistics computed, const std::string &declared, const std::string &header_type) {
    const bool by_user = declared == "min-one" || declared == "max-one";
    const bool by_header = !by_user && (header_type == "min-one" || header_type == "max-one");
    if (!by_user && !by_header) return computed;
    const std::string stated = by_user ? declared : header_type;
    if (computed.computed && ((stated == "min-one" && !computed.single) || (stated == "max-one" && !computed.uniform))) {
        computed.warning = "the mask is stated to be " + stated + (by_user ? " (-M)" : " (in the header of the masked superstring)") +
                           ", but it is not: " + (stated == "min-one" ? "a k-mer has several set occurrences" : "a k-mer of the set has an unset occurrence") +
                           "; its type is recorded as " + computed.type;
        return computed;
    }
    // Exactly known to be both: each k-mer of the set occurs once.
    const bool both = computed.computed && computed.source == "computed" && computed.single && computed.uniform;
    computed.type = both ? "min-one,max-one" : stated;
    computed.source = by_user ? "declared" : "header";
    return computed;
}

/// Build the index from the MBWT, which is consumed.
void build(index &ix, mbwt &m, int k, mask_kind mask_type, uint32_t bwt_block_bytes = 64, unsigned threads = 1) {
    ix.k = k;
    ix.size = m.bwt.size();
    ix.dollar_position = m.dollar_position;
    ix.dollar_position2 = m.dollar_position2 == mbwt::NO_ROW ? NO_ROW : m.dollar_position2;
    // A second sentinel separates the superstring from its reverse complement (doubled.h).
    ix.strands = ix.dollar_position2 == NO_ROW ? 1 : 2;
    ix.preset = ix.strands == 2 ? "fast" : "small";
    const mask_statistics statistics = classify_mask(m);
    ix.mask_type_name = statistics.type == "unknown" ? "" : statistics.type;
    ix.mask_type_source = statistics.source;
    ix.klcp = std::move(m.klcp);

    ix.bwt = quad_sequence(m.bwt.size(), bwt_block_bytes);
    std::vector<std::array<uint64_t, 4>> counted(threads, std::array<uint64_t, 4>{0, 0, 0, 0});
    // Four symbols make a byte, so the threads write to different bytes.
    parallel_chunks(m.bwt.size(), threads, parallel_chunk(4), [&](size_t begin, size_t end, unsigned thread) {
        uint64_t local[4] = {0, 0, 0, 0};
        for (size_t i = begin; i < end; ++i) {
            ix.bwt.set(i, m.bwt[i]);
            local[m.bwt[i]]++;
        }
        for (int c = 0; c < 4; ++c) counted[thread][c] += local[c];
    });
    uint64_t totals[4] = {0, 0, 0, 0};
    for (const auto &part : counted) {
        for (int c = 0; c < 4; ++c) totals[c] += part[c];
    }
    ix.bwt.finalize(threads);
    m.bwt = std::vector<byte>();
    build_phases.lap("pack_bwt");
    // The count of A includes the sentinels, which stand for the smallest suffixes.
    ix.counts[0] = 1 + (ix.dollar_position2 != NO_ROW);
    ix.counts[1] = totals[0];
    ix.counts[2] = totals[0] + totals[1];
    ix.counts[3] = totals[0] + totals[1] + totals[2];

    ix.mask_type = mask_type;
    if (mask_type == mask_kind::plain) {
        ix.mask = bit_rank(m.sa_transformed_mask.size());
        // A chunk is a whole number of blocks of the mask.
        parallel_chunks(m.sa_transformed_mask.size(), threads, parallel_chunk(bit_rank::BLOCK_BITS),
                        [&](size_t begin, size_t end, unsigned) {
                            for (size_t i = begin; i < end; ++i) {
                                if (m.sa_transformed_mask[i]) ix.mask.set(i);
                            }
                        });
        ix.mask.finalize(threads);
    } else if (mask_type == mask_kind::packed) {
        ix.mask_packed.build(m.sa_transformed_mask.size(), [&](size_t i) { return (bool)m.sa_transformed_mask[i]; });
    } else {
        ix.mask_rrr = sdsl::rrr_vector<RRR_BLOCK_SIZE>(m.sa_transformed_mask);
        ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
    }
    m.sa_transformed_mask.resize(0);
    prepare_access(ix);
    build_phases.lap("pack_mask");
}

/// The layout of an index in memory. It is chosen when the index is loaded, as the file holds the BWT and the
/// mask in a form of its own; the file recommends one (that of its preset).
struct layout {
    /// Bytes per block of the BWT (64, 128 or 256); 0 stands for the recommended one.
    uint32_t bwt_block_bytes = 0;
    /// Encoding of the mask (a mask_kind); -1 stands for the recommended one.
    int mask = -1;
};

/// The layouts of the presets. --extra-small takes the least memory: blocks of 256 bytes and the Tunstall-coded
/// mask (2.07 + 0.3 to 0.9 bits per character of the superstring). --small and --fast answer the fastest: blocks
/// of 64 bytes, i.e. one cache line per rank, and the plain mask (2.29 + 1.03 bits per character); --small
/// indexes the superstring, --fast also its reverse complement.
inline layout preset_layout(const std::string &preset) {
    layout result;
    const bool compact = preset == "extra-small";
    result.bwt_block_bytes = compact ? 256 : 64;
    result.mask = (int)(compact ? mask_kind::packed : mask_kind::plain);
    return result;
}

inline const char *mask_name(mask_kind kind) {
    return kind == mask_kind::plain ? "plain" : kind == mask_kind::packed ? "packed" : "rrr";
}

inline bool parse_mask_name(const std::string &name, mask_kind &kind) {
    if (name == "plain") {
        kind = mask_kind::plain;
    } else if (name == "packed") {
        kind = mask_kind::packed;
    } else if (name == "rrr") {
        kind = mask_kind::rrr;
    } else {
        return false;
    }
    return true;
}

/// Format 5 (index_file.h): one file with tagged sections. The sections of an index:
///   META  text, key=value per line (below)
///   BWT   encoding 1: the characters of the BWT, 2 bits each, 4 per byte (character i at bits 2 (i % 4) of
///         byte i / 4), the sentinels stored as A; not compressed, as it does not compress
///   MASK  encoding 1: the mask in the order of the suffix array, one bit per row (bit i at bit i % 8 of byte
///         i / 8); compressed
/// The kLCP array is a file of its own (<index file>.klcp) with META and
///   KLCP  encoding 1: one bit per row (bit i is set iff rows i and i + 1 share their first k - 1 characters);
///         compressed
/// The rank counts of the BWT and of the mask, the table of the first search steps and the hints of select
/// are not stored: they are computed when the file is loaded, in the layout asked for.
///
/// Keys of META of an index: kind (index), creator, k, strands (1; 2 if the reverse complement of the superstring
/// is indexed too), preset (small, fast), rows (the length of the BWT), dollar and dollar2 (the rows of the
/// sentinels), counts (the four cumulative character counts), mask_ones, mask_order (sa), layout (the recommended
/// layout in memory: bytes per block of the BWT and the mask encoding), uid (an identifier of the content, which
/// the kLCP file repeats as index_uid). Unknown keys are ignored.
static const uint32_t FORMAT_VERSION = 5;
static const uint32_t ENCODING_BWT_PACKED = 1;
static const uint32_t ENCODING_MASK_SA_BITS = 1;
/// The mask of the superstring in the order of the text, as its runs (text_mask.h); with it, the section SAMPLES
/// (encoding 1): the distance of the sampled positions, their number, and the rows of their suffixes (64-bit
/// integers).
static const uint32_t ENCODING_MASK_TEXT_RUNS = 2;
static const uint32_t ENCODING_SAMPLES = 1;
static const uint32_t ENCODING_KLCP_BITS = 1;
static const char *const INDEX_EXTENSION = ".fmsi";
static const char *const KLCP_EXTENSION = ".klcp";
/// The files of the prototype (formats 2 to 4), still read.
static const char *const OLD_INDEX_EXTENSION = ".fmsi2";

/// The mask as a bitvector, whatever its encoding.
inline sdsl::bit_vector mask_bits(const index &ix) {
    sdsl::bit_vector bits(ix.size, 0);
    for (size_t i = 0; i < ix.size; ++i) {
        if (ix.mask_get(i)) bits[i] = 1;
    }
    return bits;
}

/// What a file says about its index besides the BWT and the mask.
struct description {
    uint32_t k = 0;
    uint32_t strands = 1;
    uint64_t rows = 0;
    uint64_t dollar_position = 0;
    uint64_t dollar_position2 = NO_ROW;
    uint64_t counts[4] = {0, 0, 0, 0};
    uint64_t mask_ones = 0;
    std::string preset;
    layout recommended;
    /// The type of the mask and how it is known (mask_statistics), and the FASTA header of the masked superstring.
    std::string mask_type, mask_type_source, superstring_header;
    /// How the file stores the mask: "sa" (in the order of the suffix array) or "text" (along the superstring).
    std::string mask_order = "sa";

    /// An identifier of the content: a hash of the numbers above, which two different indexes hardly share.
    uint64_t uid() const {
        uint64_t hash = 1469598103934665603ULL;
        auto mix = [&](uint64_t value) {
            for (int byte = 0; byte < 8; ++byte) {
                hash ^= (value >> (8 * byte)) & 255;
                hash *= 1099511628211ULL;
            }
        };
        mix(k);
        mix(strands);
        mix(rows);
        mix(dollar_position);
        mix(dollar_position2);
        for (int c = 0; c < 4; ++c) mix(counts[c]);
        mix(mask_ones);
        return hash;
    }
};

inline std::string hex(uint64_t value) {
    char buffer[17];
    snprintf(buffer, sizeof(buffer), "%016llx", (unsigned long long)value);
    return buffer;
}

inline std::string creator() { return std::string("fmsi ") + VERSION; }

/// Write an index file. bwt(begin, end, buffer) gives the bytes [begin, end) of the packed BWT, from several
/// threads; mask(out) adds the section of the mask (and what belongs to it).
inline void write_index_file(const std::string &path, const description &d,
                             const std::function<void(uint64_t, uint64_t, uint8_t *)> &bwt,
                             const std::function<void(file::writer &)> &mask, unsigned threads = 1) {
    file::metadata meta;
    meta.set("kind", "index");
    meta.set("creator", creator());
    meta.set("k", d.k);
    meta.set("strands", d.strands);
    meta.set("preset", d.preset);
    meta.set("rows", d.rows);
    meta.set("dollar", d.dollar_position);
    if (d.dollar_position2 != NO_ROW) meta.set("dollar2", d.dollar_position2);
    meta.set("counts", std::to_string(d.counts[0]) + "," + std::to_string(d.counts[1]) + "," + std::to_string(d.counts[2]) +
                           "," + std::to_string(d.counts[3]));
    meta.set("mask_ones", d.mask_ones);
    meta.set("mask_type", d.mask_type.empty() ? "unknown" : d.mask_type);
    if (!d.mask_type_source.empty()) meta.set("mask_type_source", d.mask_type_source);
    if (!d.superstring_header.empty()) meta.set("superstring_header", d.superstring_header);
    meta.set("mask_order", d.mask_order);
    meta.set("layout", std::to_string(d.recommended.bwt_block_bytes) + "," + mask_name((mask_kind)d.recommended.mask));
    meta.set("uid", hex(d.uid()));
    const std::string text = meta.text();

    file::writer out(path);
    out.add("META", 1, file::REQUIRED, meta.items().size(), text.data(), text.size());
    out.add_stream("BWT", ENCODING_BWT_PACKED, file::REQUIRED, d.rows, (d.rows + 3) / 4, false, threads, bwt);
    mask(out);
    out.finish();
}

/// The section of a mask in the order of the suffix array: bits(begin, end, buffer) gives the bytes [begin, end)
/// of the bits, from several threads.
inline std::function<void(file::writer &)> mask_bits_section(uint64_t rows, unsigned threads,
                                                              std::function<void(uint64_t, uint64_t, uint8_t *)> bits) {
    return [=](file::writer &out) {
        out.add_stream("MASK", ENCODING_MASK_SA_BITS, file::REQUIRED, rows, (rows + 7) / 8, true, threads, bits);
    };
}

/// The mask of a superstring along the text, for an index file that stores it so (fmsi index --compact-file):
/// its runs (mask_runs) and its length.
struct text_mask {
    std::vector<uint8_t> runs;
    uint64_t characters = 0;
};

/// Bytes [begin, end) of a bit vector of sdsl, bit i at bit i % 8 of byte i / 8.
inline void bit_vector_bytes(const sdsl::bit_vector &bits, uint64_t begin, uint64_t end, uint8_t *buffer) {
    // The words of the vector cover its bits; a byte past them is zero.
    const uint64_t available = ((bits.size() + 63) / 64) * 8;
    const uint64_t stop = std::min(end, available);
    if (begin < stop) memcpy(buffer, reinterpret_cast<const uint8_t *>(bits.data()) + begin, stop - begin);
    if (stop < end) memset(buffer + (std::max(begin, stop) - begin), 0, end - std::max(begin, stop));
}

/// Write the index of an MBWT to a file, without building the index in memory. Returns its description.
/// With [along_text] (and the samples of the construction, mbwt::samples), the mask is stored along the superstring.
inline description save_mbwt(const mbwt &m, int k, const std::string &preset, layout recommended, const mask_statistics &statistics,
                             const std::string &path, unsigned threads = 1, const std::string &superstring_header = "",
                             const text_mask *along_text = nullptr) {
    description d;
    d.k = k;
    d.rows = m.bwt.size();
    d.dollar_position = m.dollar_position;
    d.dollar_position2 = m.dollar_position2 == mbwt::NO_ROW ? NO_ROW : m.dollar_position2;
    d.strands = d.dollar_position2 == NO_ROW ? 1 : 2;
    d.preset = preset;
    d.recommended = recommended;
    d.mask_type = statistics.type;
    d.mask_type_source = statistics.source;
    // One line of the metadata.
    d.superstring_header = superstring_header;
    std::replace(d.superstring_header.begin(), d.superstring_header.end(), '\n', ' ');
    std::vector<std::array<uint64_t, 4>> counted(threads, std::array<uint64_t, 4>{0, 0, 0, 0});
    parallel_chunks(m.bwt.size(), threads, parallel_chunk(4), [&](size_t begin, size_t end, unsigned thread) {
        uint64_t local[4] = {0, 0, 0, 0};
        for (size_t i = begin; i < end; ++i) local[m.bwt[i]]++;
        for (int c = 0; c < 4; ++c) counted[thread][c] += local[c];
    });
    uint64_t totals[4] = {0, 0, 0, 0};
    for (const auto &part : counted) {
        for (int c = 0; c < 4; ++c) totals[c] += part[c];
    }
    // The count of A includes the sentinels, which stand for the smallest suffixes.
    d.counts[0] = 1 + (d.dollar_position2 != NO_ROW);
    d.counts[1] = totals[0];
    d.counts[2] = totals[0] + totals[1];
    d.counts[3] = totals[0] + totals[1] + totals[2];
    const uint64_t *words = m.sa_transformed_mask.data();
    for (size_t w = 0; w < (m.sa_transformed_mask.size() + 63) / 64; ++w) d.mask_ones += __builtin_popcountll(words[w]);
    std::function<void(file::writer &)> mask;
    if (along_text != nullptr) {
        const uint64_t step = mbwt_sample_step, length = d.strands == 2 ? 2 * along_text->characters + 1 : along_text->characters;
        if (step == 0 || m.samples.size() != length / step) throw std::invalid_argument("the construction did not sample the suffix array");
        d.mask_order = "text";
        mask = [&m, along_text, step, threads](file::writer &out) {
            const std::vector<uint8_t> &runs = along_text->runs;
            out.add_stream("MASK", ENCODING_MASK_TEXT_RUNS, file::REQUIRED, along_text->characters, runs.size(), true, threads,
                           [&](uint64_t begin, uint64_t end, uint8_t *buffer) { memcpy(buffer, runs.data() + begin, end - begin); }, 9);
            std::vector<uint64_t> samples = {step, (uint64_t)m.samples.size()};
            samples.insert(samples.end(), m.samples.begin(), m.samples.end());
            out.add("SAMPLES", ENCODING_SAMPLES, file::REQUIRED, m.samples.size(), samples.data(), samples.size() * sizeof(uint64_t));
        };
    } else {
        mask = mask_bits_section(d.rows, threads,
                                 [&m](uint64_t begin, uint64_t end, uint8_t *buffer) { bit_vector_bytes(m.sa_transformed_mask, begin, end, buffer); });
    }
    write_index_file(
        path, d,
        [&](uint64_t begin, uint64_t end, uint8_t *buffer) {
            for (uint64_t b = begin; b < end; ++b) {
                uint8_t packed = 0;
                for (uint64_t i = 4 * b; i < std::min<uint64_t>(4 * b + 4, d.rows); ++i) packed |= m.bwt[i] << (2 * (i & 3));
                buffer[b - begin] = packed;
            }
        },
        mask, threads);
    return d;
}

inline description describe(const index &ix) {
    description d;
    d.k = ix.k;
    d.strands = ix.strands;
    d.rows = ix.size;
    d.dollar_position = ix.dollar_position;
    d.dollar_position2 = ix.dollar_position2;
    for (int c = 0; c < 4; ++c) d.counts[c] = ix.counts[c];
    d.mask_ones = ix.mask_ones;
    d.preset = ix.preset;
    d.recommended.bwt_block_bytes = ix.bwt.block_bytes();
    d.recommended.mask = (int)ix.mask_type;
    d.mask_type = ix.mask_type_name;
    d.mask_type_source = ix.mask_type_source;
    d.superstring_header = ix.superstring_header;
    return d;
}

/// Write an index that is in memory to a file; the file recommends its layout, or the given one.
inline void save(const index &ix, const std::string &path, unsigned threads = 1, const layout *recommended = nullptr) {
    description d = describe(ix);
    if (recommended != nullptr) d.recommended = *recommended;
    write_index_file(
        path, d, [&](uint64_t begin, uint64_t end, uint8_t *buffer) { ix.bwt.get_packed(begin, buffer, end - begin); },
        mask_bits_section(ix.size, threads,
                          [&ix](uint64_t begin, uint64_t end, uint8_t *buffer) {
                              for (uint64_t b = begin; b < end; ++b) {
                                  uint8_t byte = 0;
                                  for (uint64_t i = 8 * b; i < std::min<uint64_t>(8 * b + 8, ix.size); ++i) byte |= (uint8_t)ix.mask_get(i) << (i & 7);
                                  buffer[b - begin] = byte;
                              }
                          }),
        threads);
}

/// The bits of the pieces of a section, 64 at a time: sink(first bit, word, count).
template <typename sink_t>
inline void for_each_word(uint64_t offset, const uint8_t *data, size_t bytes, uint64_t total_bits, sink_t sink) {
    if (offset % 8 != 0) throw std::logic_error("pieces of a section must start at whole words");
    for (size_t b = 0; b < bytes && 8 * (offset + b) < total_bits; b += 8) {
        uint64_t word = 0;
        memcpy(&word, data + b, std::min<size_t>(8, bytes - b));
        const uint64_t first = 8 * (offset + b);
        sink(first, word, (unsigned)std::min<uint64_t>(64, total_bits - first));
    }
}

/// Load an index of the prototype (formats 2 to 4), in the layout it was stored for.
inline void load_old(index &ix, const std::string &path, unsigned threads) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::invalid_argument("couldn't open file " + path);
    char magic[sizeof(file::MAGIC)];
    in.read(magic, sizeof(magic));
    uint32_t version = 0, mask_type = 0;
    read_pod(in, version);
    if (!in || memcmp(magic, file::MAGIC, sizeof(magic)) != 0 || version < 2 || version > 4) {
        throw std::invalid_argument("file " + path + " is not an FMSI 2 index of a supported version");
    }
    read_pod(in, ix.k);
    ix.strands = 1;
    if (version >= 3) read_pod(in, ix.strands);
    read_pod(in, ix.size);
    read_pod(in, ix.dollar_position);
    ix.dollar_position2 = NO_ROW;
    if (version >= 3) read_pod(in, ix.dollar_position2);
    if (ix.strands != 1 && ix.strands != 2) throw std::invalid_argument("file " + path + " has an unsupported number of strands");
    for (int c = 0; c < 4; ++c) read_pod(in, ix.counts[c]);
    read_pod(in, mask_type);
    ix.mask_type = (mask_kind)mask_type;
    ix.preset = ix.strands == 2 ? "fast" : "small";
    if (version >= 4) {
        ix.bwt.read_packed(in, threads);
        sdsl::rrr_vector<RRR_BLOCK_SIZE> rrr;
        rrr.load(in);
        if (!in) throw std::invalid_argument("file " + path + " is truncated");
        if (ix.mask_type == mask_kind::rrr) {
            ix.mask_rrr = std::move(rrr);
            ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
        } else if (ix.mask_type == mask_kind::plain) {
            // 64 bits at a time; each thread fills whole blocks of the mask.
            ix.mask = bit_rank(ix.size);
            parallel_chunks(ix.size, threads, parallel_chunk(bit_rank::BLOCK_BITS), [&](size_t begin, size_t end, unsigned) {
                for (size_t i = begin; i < end; i += 64) {
                    const unsigned count = (unsigned)std::min<size_t>(64, end - i);
                    ix.mask.set_bits(i, rrr.get_int(i, count), count);
                }
            });
            ix.mask.finalize(threads);
        } else {
            sdsl::bit_vector bits(ix.size, 0);
            for (size_t i = 0; i < ix.size; i += 64) {
                const unsigned count = (unsigned)std::min<size_t>(64, ix.size - i);
                bits.set_int(i, rrr.get_int(i, count), count);
            }
            ix.mask_packed.build(ix.size, [&](size_t i) { return (bool)bits[i]; });
        }
    } else {
        ix.bwt.read(in);
        if (ix.mask_type == mask_kind::plain) {
            ix.mask.read(in);
        } else if (ix.mask_type == mask_kind::packed) {
            ix.mask_packed.read(in);
        } else {
            ix.mask_rrr.load(in);
            ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
        }
    }
    if (!in) throw std::invalid_argument("file " + path + " is truncated");
    prepare_access(ix);
    ix.uid = 0;
}

/// The description of an index file (format 5).
inline description read_description(const file::reader &in) {
    const file::metadata meta = in.meta();
    if (meta.get("kind") != "index") throw std::invalid_argument("file " + in.path() + " is not an index of FMSI 2");
    description d;
    d.k = (uint32_t)meta.number("k");
    d.strands = (uint32_t)meta.number("strands", 1);
    d.rows = meta.number("rows");
    d.dollar_position = meta.number("dollar");
    d.dollar_position2 = meta.number("dollar2", NO_ROW);
    d.mask_ones = meta.number("mask_ones");
    d.preset = meta.get("preset");
    d.mask_type = meta.get("mask_type", "unknown");
    d.mask_type_source = meta.get("mask_type_source");
    d.superstring_header = meta.get("superstring_header");
    std::istringstream counts(meta.get("counts"));
    std::string item;
    for (int c = 0; c < 4; ++c) {
        if (!std::getline(counts, item, ',')) throw std::invalid_argument("file " + in.path() + " has no character counts");
        d.counts[c] = std::stoull(item);
    }
    d.recommended = preset_layout(d.preset);
    const std::string recommended = meta.get("layout");
    const size_t comma = recommended.find(',');
    if (comma != std::string::npos) {
        mask_kind kind;
        const uint32_t block_bytes = (uint32_t)std::stoul(recommended.substr(0, comma));
        if (quad_sequence::valid_block_bytes(block_bytes)) d.recommended.bwt_block_bytes = block_bytes;
        if (parse_mask_name(recommended.substr(comma + 1), kind)) d.recommended.mask = (int)kind;
    }
    if (d.strands != 1 && d.strands != 2) throw std::invalid_argument("file " + in.path() + " has an unsupported number of strands");
    d.mask_order = meta.get("mask_order", "sa");
    if (d.mask_order != "sa" && d.mask_order != "text") {
        throw std::invalid_argument("file " + in.path() + " stores its mask in an order that this version of FMSI does not know");
    }
    return d;
}

/// The mask of a file that stores it along the superstring (text_mask.h), permuted into the order of the suffix
/// array by walks over the BWT, which must be loaded: into the plain mask directly, or into bits that the packed
/// mask then replaces chunk by chunk.
inline void load_text_mask(index &ix, const file::reader &in, const description &d, unsigned threads) {
    const file::section_entry &mask = in.get("MASK"), &sampled = in.get("SAMPLES");
    const uint64_t characters = mask.items;
    if (sampled.encoding != ENCODING_SAMPLES || (d.strands == 2 ? 2 * characters + 2 : characters + 1) != d.rows) {
        throw std::invalid_argument("file " + in.path() + " stores its mask in a form that this version of FMSI does not know");
    }
    // The runs and the samples are needed for the walks only: they are given back before the mask is packed.
    std::string runs = in.read_all(mask), stored = in.read_all(sampled);
    text_samples samples;
    uint64_t count = 0;
    if (stored.size() < 16) throw std::invalid_argument("file " + in.path() + " is damaged");
    memcpy(&samples.step, stored.data(), 8);
    memcpy(&count, stored.data() + 8, 8);
    if (count != (stored.size() - 16) / 8) throw std::invalid_argument("file " + in.path() + " is damaged");
    samples.rows.resize(count);
    memcpy(samples.rows.data(), stored.data() + 16, 8 * count);
    std::string().swap(stored);
    for (uint64_t row : samples.rows) {
        if (row >= d.rows) throw std::invalid_argument("file " + in.path() + " is damaged");
    }

    auto permute = [&](auto set, auto prefetch) {
        {
            const run_mask along_text(reinterpret_cast<const uint8_t *>(runs.data()), runs.size(), characters, samples.step);
            ix.bwt.with_ranks([&](const auto &ranks) {
                permute_text_mask(ranks, ix.counts, ix.dollar_position, ix.dollar_position2, ix.strands, ix.k, along_text, samples, threads, set, prefetch);
            });
        }
        std::string().swap(runs);
        std::vector<uint64_t>().swap(samples.rows);
    };
    if (ix.mask_type == mask_kind::plain) {
        ix.mask = bit_rank(d.rows);
        permute([&](uint64_t row) { ix.mask.set_atomic(row); }, [&](uint64_t row) { ix.mask.prefetch_write(row); });
        ix.mask.finalize(threads);
        return;
    }
    chunked_bits bits(d.rows);
    permute([&](uint64_t row) { bits.set_atomic(row); }, [&](uint64_t row) { bits.prefetch(row); });
    if (ix.mask_type == mask_kind::packed) {
        // The bits are passed twice; the second time they are given back as the packed mask grows.
        int passes = 0;
        ix.mask_packed.build_stream(d.rows, d.mask_ones, [&](auto sink) { bits.pass(sink, ++passes == 2); });
    } else {
        sdsl::bit_vector plain(d.rows, 0);
        uint64_t first = 0;
        bits.pass([&](uint64_t word, unsigned bit_count) {
            plain.set_int(first, word, bit_count);
            first += bit_count;
        }, true);
        ix.mask_rrr = sdsl::rrr_vector<RRR_BLOCK_SIZE>(plain);
        ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
    }
}

/// Load an index in the given layout (by default the one the file recommends). The BWT is read into its blocks
/// and the mask decoded into its encoding piece by piece, so that the peak memory is that of the loaded index;
/// the rank counts are computed on the given threads.
inline void load(index &ix, const std::string &path, unsigned threads = 1, layout wanted = layout()) {
    const uint32_t version = file::file_version(path);
    if (version >= 2 && version <= 4) {
        load_old(ix, path, threads);
        return;
    }
    file::reader in(path);
    in.check_known({"META", "BWT", "MASK", "SAMPLES"});
    const description d = read_description(in);
    ix.k = d.k;
    ix.strands = d.strands;
    ix.size = d.rows;
    ix.dollar_position = d.dollar_position;
    ix.dollar_position2 = d.dollar_position2;
    for (int c = 0; c < 4; ++c) ix.counts[c] = d.counts[c];
    ix.preset = d.preset;
    ix.uid = d.uid();
    ix.mask_type_name = d.mask_type == "unknown" ? "" : d.mask_type;
    ix.mask_type_source = d.mask_type_source;
    ix.superstring_header = d.superstring_header;
    const uint32_t block_bytes = wanted.bwt_block_bytes ? wanted.bwt_block_bytes : d.recommended.bwt_block_bytes;
    ix.mask_type = (mask_kind)(wanted.mask >= 0 ? wanted.mask : d.recommended.mask);
    if (!quad_sequence::valid_block_bytes(block_bytes)) throw std::invalid_argument("unsupported size of the blocks of the BWT");

    const file::section_entry &bwt = in.get("BWT");
    if (bwt.encoding != ENCODING_BWT_PACKED || bwt.items != d.rows || bwt.raw_bytes != (d.rows + 3) / 4) {
        throw std::invalid_argument("file " + path + " stores its BWT in a form that this version of FMSI does not know");
    }
    ix.bwt = quad_sequence(d.rows, block_bytes);
    in.read(bwt, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) { ix.bwt.set_packed(offset, data, bytes); });
    ix.bwt.finalize(threads);

    const file::section_entry &mask = in.get("MASK");
    if (mask.encoding == ENCODING_MASK_TEXT_RUNS) {
        load_text_mask(ix, in, d, threads);
    } else if (mask.encoding != ENCODING_MASK_SA_BITS || mask.items != d.rows || mask.raw_bytes != (d.rows + 7) / 8) {
        throw std::invalid_argument("file " + path + " stores its mask in a form that this version of FMSI does not know");
    } else if (ix.mask_type == mask_kind::plain) {
        ix.mask = bit_rank(d.rows);
        in.read(mask, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) {
            for_each_word(offset, data, bytes, d.rows, [&](uint64_t first, uint64_t word, unsigned count) { ix.mask.set_bits(first, word, count); });
        });
        ix.mask.finalize(threads);
    } else if (ix.mask_type == mask_kind::packed) {
        ix.mask_packed.build_stream(d.rows, d.mask_ones, [&](auto sink) {
            in.read(mask, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) {
                for_each_word(offset, data, bytes, d.rows, [&](uint64_t, uint64_t word, unsigned count) { sink(word, count); });
            });
        });
    } else {
        sdsl::bit_vector bits(d.rows, 0);
        in.read(mask, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) {
            for_each_word(offset, data, bytes, d.rows, [&](uint64_t first, uint64_t word, unsigned count) { bits.set_int(first, word, count); });
        });
        ix.mask_rrr = sdsl::rrr_vector<RRR_BLOCK_SIZE>(bits);
        ix.mask_rrr_rank = sdsl::rank_support_rrr<1, RRR_BLOCK_SIZE>(&ix.mask_rrr);
    }
    prepare_access(ix);
    if (ix.mask_ones != d.mask_ones) throw std::invalid_argument("file " + path + " is damaged: the mask does not match its description");
}

/// The kLCP array is kept in a file of its own, as it is needed for streaming only: META (kind klcp, k, rows,
/// klcp_ones, index_uid) and the bits, compressed.
inline void save_klcp_bits(const sdsl::bit_vector &klcp, uint32_t k, uint64_t uid, const std::string &path, unsigned threads = 1) {
    uint64_t ones = 0;
    for (size_t w = 0; w < (klcp.size() + 63) / 64; ++w) ones += __builtin_popcountll(klcp.data()[w]);
    file::metadata meta;
    meta.set("kind", "klcp");
    meta.set("creator", creator());
    meta.set("k", k);
    meta.set("rows", klcp.size());
    meta.set("klcp_ones", ones);
    meta.set("index_uid", hex(uid));
    const std::string text = meta.text();
    file::writer out(path);
    out.add("META", 1, file::REQUIRED, meta.items().size(), text.data(), text.size());
    out.add_stream("KLCP", ENCODING_KLCP_BITS, file::REQUIRED, klcp.size(), (klcp.size() + 7) / 8, true, threads,
                   [&](uint64_t begin, uint64_t end, uint8_t *buffer) { bit_vector_bytes(klcp, begin, end, buffer); });
    out.finish();
}

inline void save_klcp(const index &ix, const std::string &path, unsigned threads = 1) {
    save_klcp_bits(ix.klcp, ix.k, describe(ix).uid(), path, threads);
}

inline void load_klcp(index &ix, const std::string &path, unsigned threads = 1) {
    if (file::file_version(path) == 0) {
        // The kLCP file of the prototype: a bit vector of sdsl.
        if (!sdsl::load_from_file(ix.klcp, path) || ix.klcp.size() != ix.size) {
            throw std::invalid_argument("file " + path + " is not the kLCP array of the index");
        }
        return;
    }
    file::reader in(path);
    in.check_known({"META", "KLCP"});
    const file::metadata meta = in.meta();
    const file::section_entry &bits = in.get("KLCP");
    if (meta.get("kind") != "klcp" || meta.number("rows") != ix.size || bits.encoding != ENCODING_KLCP_BITS ||
        bits.items != ix.size || bits.raw_bytes != (ix.size + 7) / 8 || (ix.uid != 0 && meta.get("index_uid") != hex(ix.uid))) {
        throw std::invalid_argument("file " + path + " is not the kLCP array of the index");
    }
    ix.klcp = sdsl::bit_vector(ix.size, 0);
    in.read(bits, threads, [&](uint64_t offset, const uint8_t *data, size_t bytes) {
        for_each_word(offset, data, bytes, ix.size, [&](uint64_t first, uint64_t word, unsigned count) { ix.klcp.set_int(first, word, count); });
    });
}

/// The masked superstring of an index: its characters, in upper case where the mask is set. One walk over the
/// BWT from the end of the superstring to its start (a step per character, each a cache miss: about 0.1 s per
/// million characters of a large index).
inline std::string export_masked_superstring(const index &ix) {
    const uint64_t n = ix.strands == 2 ? (ix.size - 2) / 2 : ix.size - 1;
    std::string superstring(n, 'A');
    ix.bwt.with_ranks([&](const auto &ranks) {
        // The suffix after the last character of the superstring: the empty one, or that of the separator, which
        // is the second smallest one in an index of both strands.
        uint64_t row = ix.strands == 2 ? 1 : 0;
        for (uint64_t p = n; p-- > 0;) {
            const uint8_t c = ranks[row];
            row = ix.counts[c] + ranks.rank_inlined(row, c) - (c == 0 ? dollars_before(ix, row) : 0);
            superstring[p] = (ix.mask_get(row) ? "ACGT" : "acgt")[c];
        }
    });
    return superstring;
}

/// The file of an index given as the file itself, or as the name without its extension (a prefix, as the
/// commands of FMSI took it); an index of the prototype (<prefix>.fmsi2) is found too. Empty if there is none.
inline std::string find_index_file(const std::string &name) {
    auto is_file = [](const std::string &path) { return file::file_version(path) != 0; };
    if (is_file(name)) return name;
    if (is_file(name + INDEX_EXTENSION)) return name + INDEX_EXTENSION;
    if (is_file(name + OLD_INDEX_EXTENSION)) return name + OLD_INDEX_EXTENSION;
    return "";
}

/// The default length of the tabulated patterns, by the layout of the index: the table (8 * 4^length bytes)
/// takes at most half a bit per row of the BWT in the fast layout (blocks of 64 bytes and the plain mask) and at
/// most 1/32 bit per row in a compact one, so that it stays about one per cent of the index. The length is at
/// least 4 and at most 10, or 12 for an index of both strands in the fast layout, whose searches are the
/// shortest and gain the most from it.
inline unsigned default_prefix_length(const index &ix) {
    const bool fast = ix.bwt.block_bytes() == 64 && ix.mask_type == mask_kind::plain;
    const double budget_bits = fast ? ix.size / 2.0 : ix.size / 32.0;
    const unsigned longest = fast && ix.strands == 2 ? 12 : 10;
    unsigned length = 4;
    while (length < longest && 64.0 * ((uint64_t)1 << (2 * (length + 1))) <= budget_bits) ++length;
    return length;
}

/// Indexes of fewer rows than this get a table of 32-bit integers. A variable, so that tests reach the table of
/// 64-bit integers with small indexes.
inline uint64_t prefix_table_32bit_limit = (uint64_t)1 << 32;

/// Tabulate the SA intervals of all patterns of the given length. An entry is indexed by the
/// characters in the order of the backward search, the first one being the most significant.
void build_prefix_table(index &ix, unsigned length) {
    length = std::min(length, ix.k);
    ix.prefix_length = length;
    ix.prefix_table32.clear();
    ix.prefix_table64.clear();
    // Level by level in the integers of the table, so that the peak is the last two levels (1.25 times the table).
    auto build = [&](auto &table) {
        using value_t = typename std::decay_t<decltype(table)>::value_type;
        std::vector<value_t> current = {0, (value_t)ix.size};
        for (unsigned depth = 0; depth < length; ++depth) {
            std::vector<value_t> next(current.size() * 4);
            for (size_t key = 0; key < current.size() / 2; ++key) {
                for (byte c = 0; c < 4; ++c) {
                    size_t i = current[2 * key], j = current[2 * key + 1];
                    if (i != j) update_range(ix, i, j, c);
                    next[2 * (4 * key + c)] = (value_t)i;
                    next[2 * (4 * key + c) + 1] = (value_t)j;
                }
            }
            current = std::move(next);
        }
        table = std::move(current);
    };
    if (ix.size < prefix_table_32bit_limit) {
        build(ix.prefix_table32);
    } else {
        build(ix.prefix_table64);
    }
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
/// A k-mer is represented if the k-mer or its reverse complement has an occurrence set in the mask; if
/// the index holds both strands, the k-mer itself is searched only, and the hashes are those of the k-mers
/// as they are (a k-mer and its reverse complement have different hashes).
template <bool hashes, typename out_t>
void lookup_one_strand(const index &ix, const char *sequence, const uint32_t *offsets, size_t count, out_t *out,
                       unsigned batch_size, lookup_stats &stats);

template <bool hashes, typename out_t>
void lookup_batch(const index &ix, const char *sequence, const uint32_t *offsets, size_t count, out_t *out,
                  unsigned batch_size, lookup_stats &stats) {
    if (ix.strands == 2 && ix.k <= 32 && ix.prefix_length <= ix.k) {
        lookup_one_strand<hashes, out_t>(ix, sequence, offsets, count, out, batch_size, stats);
        return;
    }
    const unsigned k = ix.k;
    const unsigned p = ix.prefix_length;

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
    // With both strands indexed, a k-mer is searched as it is only.
    const bool one_strand = ix.strands == 2;

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
        ix.prefix_prefetch(key);
    };
    auto load = [&](slot &s) {
        if (next == count) {
            s.state = IDLE;
            --active;
            return;
        }
        s.offset = offsets[next++];
        s.kmer = sequence + s.offset;
        s.reverse = !one_strand && score < 0;
        s.second = 0;
        begin_strand(s);
    };
    auto strand_failed = [&](slot &s) {
        if (s.second || one_strand) {
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
                ix.prefix_interval(s.key, s.i, s.j);
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
                if (ix.mask_any(s.i, s.j, s.hint)) {
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

/// Lookups on an index of both strands (index2 --fast) for k <= 32, with the rank structure of the BWT
/// given (one instantiation per block size): a k-mer is searched as it is only, its symbols packed into a
/// word in the order of the search, so that a step of a slot is a shift, the ranks, and the prefetch of
/// the next block(s). Results as in lookup_batch.
template <bool hashes, typename out_t, typename ranks_t>
void lookup_one_strand(const index &ix, const ranks_t &bwt, const char *sequence, const uint32_t *offsets,
                       size_t count, out_t *out, unsigned batch_size, lookup_stats &stats) {
    const unsigned k = ix.k;
    const unsigned p = ix.prefix_length;
    const uint64_t counts[4] = {ix.counts[0], ix.counts[1], ix.counts[2], ix.counts[3]};
    const uint64_t dollar = ix.dollar_position, dollar2 = ix.dollar_position2;
    const bool mask_has_directory = ix.mask_has_directory();

    enum : uint8_t { SEARCH, START, LOCATE_MASK, MASK, IDLE };
    struct slot {
        size_t i, j;
        /// The symbols still to be searched, the next one in the lowest two bits.
        uint64_t rest;
        /// The key of the table (START), or the value returned by mask_prefetch(i) (MASK).
        uint64_t key_or_hint;
        uint32_t offset;
        /// Number of symbols still to be searched.
        uint8_t left;
        uint8_t state;
    };
    std::vector<slot> slots(std::max(1u, batch_size));
    size_t next = 0, active = slots.size();
    uint64_t found = 0, steps = 0;

    auto load = [&](slot &s) {
        if (next == count) {
            s.state = IDLE;
            --active;
            return;
        }
        s.offset = offsets[next++];
        const char *kmer = sequence + s.offset;
        // The search goes from the last symbol of the k-mer to the first one.
        uint64_t key = 0, rest = 0;
        for (unsigned t = 0; t < p; ++t) key = (key << 2) | nucleotideToInt[(uint8_t)kmer[k - 1 - t]];
        for (unsigned t = p; t < k; ++t) rest |= (uint64_t)nucleotideToInt[(uint8_t)kmer[k - 1 - t]] << (2 * (t - p));
        s.rest = rest;
        s.left = k - p;
        if (p == 0) {
            s.i = 0;
            s.j = ix.size;
            s.state = SEARCH;
        } else {
            s.key_or_hint = key;
            s.state = START;
            ix.prefix_prefetch(key);
        }
    };
    auto to_mask = [&](slot &s) {
        if (mask_has_directory) {
            s.state = LOCATE_MASK;
            ix.mask_prefetch_directory(s.i);
        } else {
            s.state = MASK;
            s.key_or_hint = ix.mask_prefetch(s.i);
        }
    };

    for (slot &s : slots) load(s);
    while (active > 0) {
        for (slot &s : slots) {
            if (__builtin_expect(s.state == SEARCH, 1)) {
                const uint8_t c = s.rest & 3;
                s.rest >>= 2;
                size_t rank_i, rank_j;
                bwt.rank_pair(s.i, s.j, c, rank_i, rank_j);
                // The sentinels are stored as A.
                const size_t dollars_i = c == 0 ? (s.i > dollar) + (s.i > dollar2) : 0;
                const size_t dollars_j = c == 0 ? (s.j > dollar) + (s.j > dollar2) : 0;
                s.i = counts[c] + rank_i - dollars_i;
                s.j = counts[c] + rank_j - dollars_j;
                ++steps;
                if (s.i == s.j) {
                    out[s.offset] = hashes ? -1 : 0;
                    load(s);
                } else if (--s.left == 0) {
                    to_mask(s);
                } else {
                    bwt.prefetch(s.i);
                    if (s.j - s.i > 1) bwt.prefetch(s.j);
                }
                continue;
            }
            switch (s.state) {
            case START:
                ix.prefix_interval(s.key_or_hint, s.i, s.j);
                if (s.i == s.j) {
                    out[s.offset] = hashes ? -1 : 0;
                    load(s);
                } else if (s.left == 0) {
                    to_mask(s);
                } else {
                    s.state = SEARCH;
                    bwt.prefetch(s.i);
                    bwt.prefetch(s.j);
                }
                break;
            case LOCATE_MASK:
                s.state = MASK;
                s.key_or_hint = ix.mask_prefetch(s.i);
                break;
            case MASK:
                if (ix.mask_any(s.i, s.j, s.key_or_hint)) {
                    out[s.offset] = hashes ? (out_t)ix.mask_rank(s.i) : (out_t)1;
                    ++found;
                } else {
                    out[s.offset] = hashes ? -1 : 0;
                }
                load(s);
                break;
            default:
                break;
            }
        }
    }
    stats.kmers += count;
    stats.found += found;
    stats.steps += steps;
    stats.searches += count;
}

template <bool hashes, typename out_t>
void lookup_one_strand(const index &ix, const char *sequence, const uint32_t *offsets, size_t count, out_t *out,
                       unsigned batch_size, lookup_stats &stats) {
    ix.bwt.with_ranks([&](const auto &ranks) {
        lookup_one_strand<hashes, out_t>(ix, ranks, sequence, offsets, count, out, batch_size, stats);
    });
}

/// A stretch of a query sequence without invalid characters, at least k characters long.
struct segment {
    uint32_t offset, length;
};

/// Look up all k-mers of the segments by streaming, [batch_size] segments at a time; the results
/// are as in lookup_batch (out must be initialized to 0, or to -1 if [hashes]).
///
/// A segment is searched on one strand and the k-mers not found there on the other one (if the index holds
/// both strands, on the strand of the query only). On a strand,
/// the k-mers (windows) are resolved from the right to the left: the match of a window is extended by one
/// character to the match of the next one, so the match can get longer than k. A longer match has a
/// subset of the occurrences of the k-mer, hence a set mask bit proves the k-mer represented, while an
/// unset one requires searching the k-mer itself. This needs no kLCP array.
///
/// If the index has the kLCP array, the last character of the match is dropped before the match is
/// extended, so the match stays a k-mer: consecutive k-mers of the query need not be consecutive in the
/// superstring and no k-mer is searched again because of an unset mask bit.
///
/// When a window cannot be found, the search stopped on a substring that does not occur in the
/// superstring, so that all windows containing it are absent on this strand. To find such substrings
/// cheaply, the unresolved windows containing the position where the search stopped are searched
/// starting from the lowest one, whose search stops the earliest.
template <bool hashes, typename out_t, bool one_strand, typename ranks_t>
void lookup_streaming(const index &ix, const ranks_t &bwt, const char *sequence, const segment *segments, size_t count,
                      out_t *out, unsigned batch_size, lookup_stats &stats) {
    const int k = ix.k;
    const int p = ix.prefix_length;
    // Values of out while a segment is processed; a non-negative hash or 1 means found.
    const out_t UNRESOLVED = hashes ? -1 : 0;
    const out_t ABSENT[2] = {(out_t)(hashes ? -2 : 2), (out_t)(hashes ? -3 : 3)};
    const int NO_BLOCK = -1;
    const bool klcp = ix.has_klcp();

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
    // Go from the interval of P to that of c P, with the rank structure given.
    auto extend = [&](size_t &i, size_t &j, byte c) FMSI2_INLINE_LAMBDA {
        size_t rank_i, rank_j;
        bwt.rank_pair_near(i, j, c, rank_i, rank_j);
        size_t dollars_i = c == 0 ? dollars_before(ix, i) : 0, dollars_j = c == 0 ? dollars_before(ix, j) : 0;
        i = ix.counts[c] + rank_i - dollars_i;
        j = ix.counts[c] + rank_j - dollars_j;
    };
    // The same for the interval [i, i + 1) of a single row: the row is followed if its symbol is c (and
    // not a sentinel, stored as A).
    auto extend_single = [&](size_t &i, size_t &j, byte c) FMSI2_INLINE_LAMBDA {
        uint8_t symbol;
        const size_t rank_i = bwt.rank_symbol(i, c, symbol);
        const bool follows = symbol == c && !(c == 0 && (i == ix.dollar_position || i == ix.dollar_position2));
        const size_t dollars_i = c == 0 ? dollars_before(ix, i) : 0;
        i = ix.counts[c] + rank_i - dollars_i;
        j = i + follows;
    };
    const bool plain_mask = ix.mask_type == mask_kind::plain;
    std::vector<slot> slots(std::max(1u, batch_size));
    size_t next = 0, active = slots.size();
    int score = 0;
    uint64_t found = 0, steps = 0, kmers = 0, searches = 0, rechecks = 0;

    // Character of the strand at position v, and the result of its window y.
    auto symbol = [&](const slot &s, int v) FMSI2_INLINE_LAMBDA -> byte {
        return !one_strand && s.reverse ? 3 - nucleotideToInt[(uint8_t)s.sequence[s.length - 1 - v]]
                         : nucleotideToInt[(uint8_t)s.sequence[v]];
    };
    auto result = [&](const slot &s, int y) FMSI2_INLINE_LAMBDA -> out_t & { return s.out[!one_strand && s.reverse ? s.windows - 1 - y : y]; };
    auto unresolved = [&](const slot &s, int y) FMSI2_INLINE_LAMBDA {
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
        ix.prefix_prefetch(key);
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
        if (!s.second && !one_strand) {
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
        s.reverse = !one_strand && score < 0;
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
    auto to_check = [&](slot &s) FMSI2_INLINE_LAMBDA {
        s.state = CHECK;
        if (plain_mask) {
            ix.mask.prefetch(s.i);
        } else {
            s.hint = ix.mask_prefetch(s.i);
        }
        bwt.prefetch(s.i);
        if (s.j - s.i > 1) {
            bwt.prefetch(s.j);
            if (klcp) ix.klcp_prefetch(s.i, s.j);
        } else if (klcp) {
            // The kLCP bits of rows i - 1 and i, in one word unless i starts one.
            ix.klcp_prefetch(s.i, s.i);
            if ((s.i & 63) == 0) ix.klcp_prefetch(s.j, s.j);
        }
    };
    auto after_search_step = [&](slot &s) FMSI2_INLINE_LAMBDA {
        if (s.i == s.j) {
            search_failed(s);
        } else if (s.t == k) {
            s.exact = 1;
            to_check(s);
        } else {
            s.state = SEARCH;
            bwt.prefetch(s.i);
            bwt.prefetch(s.j);
        }
    };
    // The interval is the one of a string starting with window x: resolve the window and go to the next one.
    auto check = [&](slot &s) FMSI2_INLINE_LAMBDA {
        if (plain_mask && s.j == s.i + 1 ? ix.mask[s.i] : ix.mask_any(s.i, s.j, s.hint)) {
            // No bit is set between the start of the interval and its first set bit.
            result(s, s.x) = hashes ? (out_t)ix.mask_rank(s.i) : (out_t)1;
            ++found;
            if (!one_strand) ++s.hits[s.reverse];
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
        // With the kLCP array, drop the last character of the k-mer first, so that the match stays a k-mer. Most
        // matches are a single row that this leaves as it is.
        if (s.j == s.i + 1 && (!klcp || ix.klcp_single(s.i))) {
            extend_single(s.i, s.j, symbol(s, s.x));
        } else {
            if (klcp) ix.klcp_widen(s.i, s.j);
            extend(s.i, s.j, symbol(s, s.x));
        }
        ++steps;
        if (s.i != s.j) {
            s.exact = klcp;
            to_check(s);
        } else if (klcp) {
            // The window does not occur although its last k - 1 characters do: searching it would tell nothing more.
            s.t = k;
            search_failed(s);
        } else {
            begin_search(s, s.x);
        }
    };

    for (slot &s : slots) load(s);
    while (active > 0) {
        for (slot &s : slots) {
            switch (s.state) {
            case START:
                ix.prefix_interval(s.key, s.i, s.j);
                s.t = p;
                after_search_step(s);
                break;
            case SEARCH:
                if (s.j == s.i + 1) {
                    extend_single(s.i, s.j, symbol(s, s.x + k - 1 - s.t));
                } else {
                    extend(s.i, s.j, symbol(s, s.x + k - 1 - s.t));
                }
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

/// Streamed lookups (see the instantiations below), for the rank structure of the block size of the index and
/// with the strand logic compiled out if the index holds both strands.
template <bool hashes, typename out_t>
void lookup_streaming(const index &ix, const char *sequence, const segment *segments, size_t count, out_t *out,
                      unsigned batch_size, lookup_stats &stats) {
    ix.bwt.with_ranks([&](const auto &ranks) {
        if (ix.strands == 2) {
            lookup_streaming<hashes, out_t, true>(ix, ranks, sequence, segments, count, out, batch_size, stats);
        } else {
            lookup_streaming<hashes, out_t, false>(ix, ranks, sequence, segments, count, out, batch_size, stats);
        }
    });
}

/// The first symbol of the suffix of a row (not the row 0 of the empty suffix).
inline uint8_t first_symbol(const index &ix, size_t row) {
    return (row >= ix.counts[1]) + (row >= ix.counts[2]) + (row >= ix.counts[3]);
}

/// For a row whose suffix starts with c: the number, among the occurrences of c stored in the BWT, of
/// the occurrence at the row of the suffix shortened by this c. Selecting it inverts update_range.
inline size_t next_row_rank(const index &ix, size_t row, uint8_t c) {
    size_t r = row - ix.counts[c];
    if (c != 0) return r;
    // The sentinels are stored as A in the BWT without being occurrences of A.
    r += r >= ix.dollar_rank;
    return r + (r >= ix.dollar_rank2);
}

/// Access: the k-mer with the given hash, i.e., the k-mer at the set occurrence number [hash] (from 0)
/// of the mask in the order of the suffix array; k characters are written to out. The hash must be less
/// than mask_ones. With one set occurrence per k-mer (a mask with the minimum number of ones), this
/// inverts the hashes of the lookups: the k-mer is returned as it occurs in the superstring, so it is
/// the looked-up k-mer or its reverse complement.
///
/// The row of the occurrence is found by a select on the mask, and each character by a select on the BWT.
inline void access(const index &ix, uint64_t hash, char *out) {
    size_t row = ix.mask_select(hash);
    for (uint32_t t = 0; t < ix.k; ++t) {
        uint8_t c = first_symbol(ix, row);
        out[t] = "ACGT"[c];
        if (t + 1 < ix.k) row = ix.bwt.select(next_row_rank(ix, row, c), c);
    }
}

/// Access for [count] hashes, [batch_size] of them at a time: the k-mer of hashes[i] is written to
/// out + i * k. Each select is split into its steps, and the other k-mers of the batch are advanced
/// while the memory it prefetched arrives: the hints of the select (if built), then the blocks.
inline void access_batch(const index &ix, const uint64_t *hashes, size_t count, char *out, unsigned batch_size) {
    const uint32_t k = ix.k;
    const bool plain = ix.mask_type == mask_kind::plain;
    const bool hints = ix.bwt.has_select_hints();
    // FETCH: the slot needs a hash; MASK: a select on the mask was prepared; HINT, SYMBOL: the hints
    // or the blocks of a select on the BWT were prefetched.
    enum : uint8_t { FETCH, MASK, HINT, SYMBOL, DONE };
    struct slot {
        select_state select;
        size_t row, rank;
        char *out;
        uint32_t t;
        uint8_t c, state;
    };
    std::vector<slot> slots(std::max(1u, batch_size));
    for (slot &s : slots) s.state = FETCH;
    size_t next = 0, active = slots.size();

    while (active > 0) {
        for (slot &s : slots) {
            if (s.state == DONE) continue;
            if (s.state == HINT) {
                s.select = ix.bwt.select_prepare(s.rank, s.c);
                s.state = SYMBOL;
                continue;
            }
            if (s.state == MASK) {
                s.row = ix.mask.select_finish(s.select);
            } else if (s.state == SYMBOL) {
                s.row = ix.bwt.select_finish(s.select, s.c);
            }
            // Go on until the next select is under way or no hash is left.
            for (;;) {
                if (s.state == FETCH) {
                    if (next == count) {
                        s.state = DONE;
                        --active;
                        break;
                    }
                    s.out = out + next * k;
                    s.t = 0;
                    uint64_t hash = hashes[next++];
                    if (plain) {
                        s.select = ix.mask.select_prepare(hash);
                        s.state = MASK;
                        break;
                    }
                    s.row = ix.mask_select(hash);
                }
                s.c = first_symbol(ix, s.row);
                s.out[s.t] = "ACGT"[s.c];
                if (++s.t == k) {
                    s.state = FETCH;
                    continue;
                }
                s.rank = next_row_rank(ix, s.row, s.c);
                if (hints) {
                    ix.bwt.select_prefetch_hints(s.rank, s.c);
                    s.state = HINT;
                } else {
                    s.select = ix.bwt.select_prepare(s.rank, s.c);
                    s.state = SYMBOL;
                }
                break;
            }
        }
    }
}

} // namespace fmsi2
