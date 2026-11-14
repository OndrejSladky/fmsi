#pragma once

// Rank structures of the FMSI 2 prototype.
//
// Both structures store the data and the rank counts of a block in the same
// 64-byte cache line, so that a rank query touches a single line (plus a small
// superblock array that stays in the cache). This follows the one-line rank of
// SBWT-squeeze (Diseth et al., 2026) and QuadRank (Groot Koerkamp, 2026).

#include "parallel.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <vector>

// FMSI2_SCALAR selects the portable code also on ARM (for testing).
#if defined(__aarch64__) && !defined(FMSI2_SCALAR)
#include <arm_neon.h>
#define FMSI2_NEON 1
#endif

#ifndef FMSI2_HINT_STEP
#define FMSI2_HINT_STEP 128
#endif

namespace fmsi2 {

/// A zero-initialized heap array aligned to 256 bytes (the largest block, and a multiple of a cache line).
template <typename T>
class aligned_array {
    T *data_ = nullptr;
    size_t size_ = 0;

  public:
    aligned_array() = default;
    explicit aligned_array(size_t size) { resize(size); }
    aligned_array(const aligned_array &) = delete;
    aligned_array &operator=(const aligned_array &) = delete;
    aligned_array(aligned_array &&other) noexcept : data_(other.data_), size_(other.size_) {
        other.data_ = nullptr;
        other.size_ = 0;
    }
    aligned_array &operator=(aligned_array &&other) noexcept {
        if (this != &other) {
            free(data_);
            data_ = other.data_;
            size_ = other.size_;
            other.data_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }
    ~aligned_array() { free(data_); }

    void resize(size_t size) {
        free(data_);
        data_ = nullptr;
        size_ = size;
        if (size == 0) return;
        void *memory = nullptr;
        if (posix_memalign(&memory, 256, size * sizeof(T)) != 0) throw std::bad_alloc();
        memset(memory, 0, size * sizeof(T));
        data_ = static_cast<T *>(memory);
    }

    T &operator[](size_t i) { return data_[i]; }
    const T &operator[](size_t i) const { return data_[i]; }
    T *data() { return data_; }
    const T *data() const { return data_; }
    size_t size() const { return size_; }
};

template <typename T>
void write_pod(std::ostream &out, const T &value) {
    out.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template <typename T>
void read_pod(std::istream &in, T &value) {
    in.read(reinterpret_cast<char *>(&value), sizeof(T));
}

/// The low [bits] bits set, for bits in (-inf, +inf) clamped to [0, 64].
inline uint64_t low_bits(int bits) {
    if (bits >= 64) return ~0ULL;
    if (bits <= 0) return 0ULL;
    return (1ULL << bits) - 1;
}

/// The first x in [lo, hi) for which above(x) holds, or hi if there is none; above must be monotone
/// (false up to some point, true from it on). The search gallops from the guess, so it takes few
/// evaluations when the guess is near.
template <typename predicate>
inline size_t first_above(size_t lo, size_t hi, size_t guess, predicate above) {
    if (lo >= hi) return hi;
    guess = std::min(std::max(guess, lo), hi - 1);
    size_t step = 1;
    if (above(guess)) {
        hi = guess;
        while (lo < hi) {
            size_t probe = hi - lo > step ? hi - step : lo;
            if (!above(probe)) {
                lo = probe + 1;
                break;
            }
            hi = probe;
            step *= 2;
        }
    } else {
        lo = guess + 1;
        while (lo < hi) {
            size_t probe = hi - lo > step ? lo + step - 1 : hi - 1;
            if (above(probe)) {
                hi = probe;
                break;
            }
            lo = probe + 1;
            step *= 2;
        }
    }
    while (lo < hi) {
        size_t middle = lo + (hi - lo) / 2;
        if (above(middle)) {
            hi = middle;
        } else {
            lo = middle + 1;
        }
    }
    return lo;
}

/// Position of the set bit number q (from 0) of the word, which has more than q set bits.
inline int select_in_word(uint64_t word, size_t q) {
    for (; q > 0; --q) word &= word - 1;
    return __builtin_ctzll(word);
}

/// A select query between its two steps: select_prepare finds the superblock, whose counts stay in
/// the cache, and prefetches the block where the occurrence is expected; select_finish reads the blocks.
/// In between, other queries can be advanced while the memory arrives.
struct select_state {
    /// The rank, the range of blocks to search, and the block to read first.
    size_t rank, first, last, guess;
};

/// Rank over a sequence of 2-bit symbols in blocks of BYTES bytes: 64 (2.29 bits per symbol), 128 (2.13)
/// or 256 (2.06).
///
/// A block holds 4 * (BYTES - 8) symbols and, in its middle, the numbers of occurrences of the four
/// symbols before the middle (16 bits each, relative to the superblock). A rank counts from the middle
/// to the position, so it scans at most half of the block, and touches only the cache lines between the
/// position and the counts: one line of a 64-byte block, 1.5 lines of a 128-byte block on average and 2
/// of a 256-byte one.
///
/// Layout: BYTES / 2 - 8 bytes of symbols (the left half), 8 bytes of counts, BYTES / 2 bytes of symbols
/// (the right half). Symbol t of the block is in byte t / 4 (left half) or 8 + t / 4 (right half), at
/// bits 2 * (t % 4).
template <size_t BYTES>
class quad_rank {
  public:
    static constexpr size_t BLOCK_BYTES = BYTES;
    static constexpr size_t BLOCK_SYMBOLS = 4 * (BYTES - 8);
    /// Number of symbols of the left half, and the byte of the counts.
    static constexpr size_t LEFT_SYMBOLS = 4 * (BYTES / 2 - 8);
    static constexpr size_t COUNTS_BYTE = BYTES / 2 - 8;
    /// A superblock has less than 2^16 symbols, so that the relative counts fit into 16 bits.
    static constexpr size_t SUPER_BLOCKS = BYTES <= 64 ? 256 : BYTES <= 128 ? 128 : 64;

  private:
    struct alignas(64) block_t {
        uint8_t bytes[BYTES];
    };
    static_assert(sizeof(block_t) == BYTES && BYTES % 64 == 0, "blocks must be whole cache lines");
    static_assert(SUPER_BLOCKS * BLOCK_SYMBOLS < 65536 && BYTES <= 256, "counts must fit");

    struct super_block {
        uint64_t counts[4];
    };

    /// The offsets 0, ..., BYTES / 2 - 1 of the bytes of a half.
    struct offset_table {
        alignas(64) uint8_t offsets[BYTES / 2];
        offset_table() {
            for (size_t byte = 0; byte < BYTES / 2; ++byte) offsets[byte] = (uint8_t)byte;
        }
    };
    static inline const offset_table HALF{};

    size_t size_ = 0;
    aligned_array<block_t> blocks_;
    std::vector<super_block> supers_;

    static inline size_t byte_of(size_t r) { return r / 4 + (r >= LEFT_SYMBOLS ? 8 : 0); }

    static inline uint16_t count_before_middle(const block_t &block, uint8_t c) {
        uint16_t count;
        memcpy(&count, block.bytes + COUNTS_BYTE + 2 * c, 2);
        return count;
    }

    /// The bits 2t of the byte are set for its symbols equal to c.
    static inline uint8_t matches(uint8_t byte, uint8_t c) {
        uint8_t x = byte ^ (uint8_t)(0x55 * (c ^ 3));
        return x & (x >> 1) & 0x55;
    }

    /// Number of occurrences of c among the first r symbols of the block, relative to its middle
    /// (negative in the left half). Branch-free on purpose: the position in the block is unpredictable.
    static inline int64_t count_from_middle(const block_t &block, size_t r, uint8_t c) {
        const bool right = r >= LEFT_SYMBOLS;
        // The half is scanned as BYTES / 2 bytes starting at [base]; in the left half, these include the
        // counts, which are never selected. The selected bytes are [first, first + selected) of them: from
        // the byte with symbol r to the middle (left), or from the middle to the byte with symbol r (right).
        const uint8_t *base = block.bytes + (right ? BYTES / 2 : 0);
        const size_t first = right ? 0 : r / 4;
        const size_t selected = right ? (r - LEFT_SYMBOLS) / 4 : COUNTS_BYTE - r / 4;
        int64_t sum;
#ifdef FMSI2_NEON
        const uint8x16_t pattern = vdupq_n_u8(0x55 * (c ^ 3)), pairs = vdupq_n_u8(0x55);
        const uint8x16_t first_byte = vdupq_n_u8(first), selected_bytes = vdupq_n_u8(selected);
        uint8x16_t sums = vdupq_n_u8(0);
        for (size_t v = 0; v < BYTES / 32; ++v) {
            uint8x16_t x = veorq_u8(vld1q_u8(base + 16 * v), pattern);
            uint8x16_t ones = vcntq_u8(vandq_u8(vandq_u8(x, vshrq_n_u8(x, 1)), pairs));
            uint8x16_t in_range = vcltq_u8(vsubq_u8(vld1q_u8(HALF.offsets + 16 * v), first_byte), selected_bytes);
            sums = vaddq_u8(sums, vandq_u8(ones, in_range));
        }
        sum = vaddlvq_u8(sums);
#else
        sum = 0;
        for (size_t byte = first; byte < first + selected; ++byte) sum += __builtin_popcount(matches(base[byte], c));
#endif
        // The symbols before position r in its byte: counted above in the left half, missing in the right one.
        uint8_t partial = matches(block.bytes[byte_of(r)], c) & ((1u << (2 * (r % 4))) - 1);
        int64_t before = (partial & 1) + ((partial >> 2) & 1) + ((partial >> 4) & 1);
        return right ? sum + before : before - sum;
    }

  public:
    quad_rank() = default;
    explicit quad_rank(size_t size) : size_(size), blocks_(size / BLOCK_SYMBOLS + 1) {}

    size_t size() const { return size_; }

    size_t bytes() const { return blocks_.size() * sizeof(block_t) + supers_.size() * sizeof(super_block); }

    /// Set symbol i; valid only before finalize() and only once per position.
    void set(size_t i, uint8_t c) {
        size_t r = i % BLOCK_SYMBOLS;
        blocks_[i / BLOCK_SYMBOLS].bytes[byte_of(r)] |= c << (2 * (r % 4));
    }

    /// Compute the rank counts; to be called once after all symbols are set.
    void finalize(unsigned threads = 1) {
        // The counts of a block are relative to its superblock, so the superblocks are counted independently
        // and their totals are summed up afterwards.
        size_t count = (blocks_.size() + SUPER_BLOCKS - 1) / SUPER_BLOCKS;
        supers_.assign(blocks_.size() / SUPER_BLOCKS + 1, super_block{{0, 0, 0, 0}});
        std::vector<super_block> sums(count, super_block{{0, 0, 0, 0}});
        size_t chunk = std::max<size_t>(1, parallel_chunk_items / (SUPER_BLOCKS * BLOCK_SYMBOLS));
        parallel_chunks(count, threads, chunk, [&](size_t begin, size_t end, unsigned) {
            for (size_t s = begin; s < end; ++s) {
                uint64_t totals[4] = {0, 0, 0, 0};
                size_t last = std::min(blocks_.size(), (s + 1) * SUPER_BLOCKS);
                for (size_t b = s * SUPER_BLOCKS; b < last; ++b) {
                    for (size_t r = 0; r < BLOCK_SYMBOLS; ++r) {
                        if (r == LEFT_SYMBOLS) {
                            for (int c = 0; c < 4; ++c) {
                                uint16_t relative = (uint16_t)totals[c];
                                memcpy(blocks_[b].bytes + COUNTS_BYTE + 2 * c, &relative, 2);
                            }
                        }
                        totals[(blocks_[b].bytes[byte_of(r)] >> (2 * (r % 4))) & 3]++;
                    }
                }
                for (int c = 0; c < 4; ++c) sums[s].counts[c] = totals[c];
            }
        });
        uint64_t totals[4] = {0, 0, 0, 0};
        for (size_t s = 0; s < count; ++s) {
            for (int c = 0; c < 4; ++c) {
                supers_[s].counts[c] = totals[c];
                totals[c] += sums[s].counts[c];
            }
        }
    }

    inline uint8_t operator[](size_t i) const {
        size_t r = i % BLOCK_SYMBOLS;
        return (blocks_[i / BLOCK_SYMBOLS].bytes[byte_of(r)] >> (2 * (r % 4))) & 3;
    }

    /// Prefetch the cache lines read by rank(i, c): those from the position to the counts.
    inline void prefetch(size_t i) const {
        const uint8_t *block = blocks_[i / BLOCK_SYMBOLS].bytes;
        size_t line = byte_of(i % BLOCK_SYMBOLS) / 64, counts_line = COUNTS_BYTE / 64;
        __builtin_prefetch(block + 64 * counts_line);
        if (BYTES > 64) __builtin_prefetch(block + 64 * line);
        if (BYTES > 128) __builtin_prefetch(block + 64 * ((line + counts_line + 1) / 2));
    }

    /// Number of occurrences of c among the first i symbols.
    inline size_t rank(size_t i, uint8_t c) const {
        size_t b = i / BLOCK_SYMBOLS;
        const block_t &block = blocks_[b];
        return supers_[b / SUPER_BLOCKS].counts[c] + count_before_middle(block, c) +
               count_from_middle(block, i % BLOCK_SYMBOLS, c);
    }

    /// Ranks of c at the positions i <= j; cheaper than two rank queries for neighbouring positions.
    inline void rank_pair(size_t i, size_t j, uint8_t c, size_t &rank_i, size_t &rank_j) const {
        rank_i = rank(i, c);
        rank_j = j == i + 1 ? rank_i + ((*this)[i] == c) : rank(j, c);
    }

  private:
    static constexpr size_t NOT_FOUND = ~(size_t)0;

    /// Offset of the occurrence number q (from 0) of c among the symbols of [count] bytes (a multiple of 8),
    /// or NOT_FOUND, with q decreased by the number of occurrences.
    static inline size_t select_in_half(const uint8_t *bytes, size_t count, uint8_t c, size_t &q) {
        const uint64_t pairs = 0x5555555555555555ULL, pattern = pairs * (uint64_t)(c ^ 3);
        for (size_t byte = 0; byte < count; byte += 8) {
            uint64_t word;
            memcpy(&word, bytes + byte, 8);
            uint64_t x = word ^ pattern;
            uint64_t found = x & (x >> 1) & pairs;
            size_t occurrences = __builtin_popcountll(found);
            if (q < occurrences) return 4 * byte + select_in_word(found, q) / 2;
            q -= occurrences;
        }
        return NOT_FOUND;
    }

    /// Select hints, built on demand and not stored: for each symbol, the block that select reads for every
    /// HINT_STEP-th occurrence (the last block with at most that many occurrences before its middle). They
    /// take 1/16 bit per symbol and bound the search to a few blocks (4 on average for blocks of 64 bytes),
    /// which are prefetched together.
    static constexpr size_t HINT_STEP = FMSI2_HINT_STEP;
    /// Ranges of at most this many blocks are prefetched whole.
    static constexpr size_t PREFETCHED_RANGE = 6;
    std::vector<uint32_t> hints_[4];

    /// Number of occurrences of c before the middle of block b.
    inline size_t before_middle(size_t b, uint8_t c) const {
        return supers_[b / SUPER_BLOCKS].counts[c] + count_before_middle(blocks_[b], c);
    }

  public:
    void build_select_hints() {
        size_t next[4] = {0, 0, 0, 0};
        for (int c = 0; c < 4; ++c) hints_[c].clear();
        for (size_t b = 0; b < blocks_.size(); ++b) {
            for (int c = 0; c < 4; ++c) {
                size_t middle = before_middle(b, c);
                for (; next[c] < middle; next[c] += HINT_STEP) hints_[c].push_back(b == 0 ? 0 : b - 1);
            }
        }
        // The occurrences after the last middle, and an entry to interpolate towards.
        for (int c = 0; c < 4; ++c) hints_[c].insert(hints_[c].end(), 2, (uint32_t)(blocks_.size() - 1));
    }

    bool has_select_hints() const { return !hints_[0].empty(); }

    size_t select_hint_bytes() const {
        size_t entries = 0;
        for (int c = 0; c < 4; ++c) entries += hints_[c].size();
        return entries * sizeof(uint32_t);
    }

    /// Before the first step of select(r, c): prefetch the hints it reads, if there are any.
    inline void select_prefetch_hints(size_t r, uint8_t c) const {
        if (!hints_[c].empty()) __builtin_prefetch(hints_[c].data() + std::min(r / HINT_STEP, hints_[c].size() - 2));
    }

    /// First step of select(r, c).
    inline select_state select_prepare(size_t r, uint8_t c) const {
        select_state state;
        state.rank = r;
        const std::vector<uint32_t> &hints = hints_[c];
        if (!hints.empty()) {
            // The block is between those of the neighbouring hints.
            size_t j = std::min(r / HINT_STEP, hints.size() - 2);
            state.first = hints[j];
            state.last = std::min(blocks_.size(), (size_t)hints[j + 1] + 2);
            state.guess = state.first + (r - j * HINT_STEP) * (hints[j + 1] - hints[j] + 1) / HINT_STEP;
        } else {
            // The last superblock with at most r occurrences before it; the guesses assume them evenly spread.
            const size_t count = (blocks_.size() + SUPER_BLOCKS - 1) / SUPER_BLOCKS;
            const uint64_t before_last = supers_[count - 1].counts[c];
            size_t guess = before_last > 0 ? (size_t)((unsigned __int128)r * (count - 1) / before_last) : 0;
            size_t s = first_above(1, count, guess + 1, [&](size_t x) { return supers_[x].counts[c] > r; }) - 1;
            state.first = s * SUPER_BLOCKS;
            state.last = std::min(blocks_.size(), state.first + SUPER_BLOCKS);
            size_t blocks = state.last - state.first, rest = r - supers_[s].counts[c];
            size_t in_super =
                s + 1 < count ? supers_[s + 1].counts[c] - supers_[s].counts[c] : blocks * BLOCK_SYMBOLS / 4;
            state.guess = state.first + rest * blocks / std::max<size_t>(1, in_super);
        }
        state.guess = std::min(state.guess, state.last - 1);
        // The occurrence is after the middle of a block or before the middle of the next one. With the
        // hints, the whole range is prefetched if it is short; else the blocks around the guessed one, as
        // the guess is often off by one.
        size_t from = state.guess > 0 ? state.guess - 1 : 0, to = state.guess + 2;
        if (!hints.empty() && state.last - state.first <= PREFETCHED_RANGE) {
            from = state.first;
            to = state.last - 1;
        }
        const uint8_t *end = blocks_[std::min(to, blocks_.size() - 1)].bytes + COUNTS_BYTE;
        for (const uint8_t *line = blocks_[from].bytes + COUNTS_BYTE; line <= end; line += 64) __builtin_prefetch(line);
        return state;
    }

    /// Second step of select(r, c).
    inline size_t select_finish(const select_state &state, uint8_t c) const {
        // The first block of the range with more than r occurrences before its middle.
        size_t above = first_above(state.first, state.last, state.guess,
                                   [&](size_t b) { return before_middle(b, c) > state.rank; });
        size_t q = state.rank;
        if (above > state.first) {
            q -= before_middle(above - 1, c);
            size_t symbol = select_in_half(blocks_[above - 1].bytes + BYTES / 2, BYTES / 2, c, q);
            if (symbol != NOT_FOUND) return (above - 1) * BLOCK_SYMBOLS + LEFT_SYMBOLS + symbol;
        } else {
            // Before the first middle of the range: the range starts at a superblock or at the first block.
            q -= supers_[above / SUPER_BLOCKS].counts[c];
        }
        return above * BLOCK_SYMBOLS + select_in_half(blocks_[above].bytes, COUNTS_BYTE, c, q);
    }

    /// Position of the occurrence number r (from 0) of c; there must be more than r occurrences.
    inline size_t select(size_t r, uint8_t c) const { return select_finish(select_prepare(r, c), c); }

    void write(std::ostream &out) const {
        write_pod(out, (uint64_t)size_);
        write_pod(out, (uint64_t)blocks_.size());
        out.write(reinterpret_cast<const char *>(blocks_.data()), blocks_.size() * sizeof(block_t));
        write_pod(out, (uint64_t)supers_.size());
        out.write(reinterpret_cast<const char *>(supers_.data()), supers_.size() * sizeof(super_block));
    }

    void read(std::istream &in) {
        uint64_t size, block_count, super_count;
        read_pod(in, size);
        read_pod(in, block_count);
        size_ = size;
        blocks_.resize(block_count);
        in.read(reinterpret_cast<char *>(blocks_.data()), block_count * sizeof(block_t));
        read_pod(in, super_count);
        supers_.resize(super_count);
        in.read(reinterpret_cast<char *>(supers_.data()), super_count * sizeof(super_block));
    }
};

/// A sequence of 2-bit symbols with rank, in blocks of 64, 128 or 256 bytes, chosen when it is built.
/// Larger blocks have fewer counts per symbol (2.29, 2.13 and 2.06 bits per symbol), but a rank reads
/// more cache lines and counts over more symbols.
class quad_sequence {
    uint32_t block_bytes_ = 64;
    quad_rank<64> small_;
    quad_rank<128> medium_;
    quad_rank<256> large_;

    template <typename function>
    inline auto visit(function f) const {
        switch (block_bytes_) {
        case 64:
            return f(small_);
        case 128:
            return f(medium_);
        default:
            return f(large_);
        }
    }

    template <typename function>
    inline auto visit(function f) {
        switch (block_bytes_) {
        case 64:
            return f(small_);
        case 128:
            return f(medium_);
        default:
            return f(large_);
        }
    }

  public:
    quad_sequence() = default;

    /// Call f with the rank structure of the chosen block size, for loops instantiated for each of them.
    template <typename function>
    inline auto with_ranks(function f) const {
        return visit(f);
    }

    quad_sequence(size_t size, uint32_t block_bytes) : block_bytes_(block_bytes) {
        visit([&](auto &ranks) { ranks = std::decay_t<decltype(ranks)>(size); });
    }

    static bool valid_block_bytes(uint32_t bytes) { return bytes == 64 || bytes == 128 || bytes == 256; }

    uint32_t block_bytes() const { return block_bytes_; }

    size_t bytes() const {
        return visit([](const auto &ranks) { return ranks.bytes(); });
    }

    void set(size_t i, uint8_t c) {
        visit([&](auto &ranks) { ranks.set(i, c); });
    }

    void finalize(unsigned threads = 1) {
        visit([&](auto &ranks) { ranks.finalize(threads); });
    }

    inline void prefetch(size_t i) const {
        visit([&](const auto &ranks) { ranks.prefetch(i); });
    }

    inline void rank_pair(size_t i, size_t j, uint8_t c, size_t &rank_i, size_t &rank_j) const {
        visit([&](const auto &ranks) { ranks.rank_pair(i, j, c, rank_i, rank_j); });
    }

    inline size_t rank(size_t i, uint8_t c) const {
        return visit([&](const auto &ranks) { return ranks.rank(i, c); });
    }

    void build_select_hints() {
        visit([&](auto &ranks) { ranks.build_select_hints(); });
    }

    bool has_select_hints() const {
        return visit([](const auto &ranks) { return ranks.has_select_hints(); });
    }

    size_t select_hint_bytes() const {
        return visit([](const auto &ranks) { return ranks.select_hint_bytes(); });
    }

    inline void select_prefetch_hints(size_t r, uint8_t c) const {
        visit([&](const auto &ranks) { ranks.select_prefetch_hints(r, c); });
    }

    inline select_state select_prepare(size_t r, uint8_t c) const {
        return visit([&](const auto &ranks) { return ranks.select_prepare(r, c); });
    }

    inline size_t select_finish(const select_state &state, uint8_t c) const {
        return visit([&](const auto &ranks) { return ranks.select_finish(state, c); });
    }

    inline size_t select(size_t r, uint8_t c) const {
        return visit([&](const auto &ranks) { return ranks.select(r, c); });
    }

    void write(std::ostream &out) const {
        write_pod(out, block_bytes_);
        visit([&](const auto &ranks) { ranks.write(out); });
    }

    void read(std::istream &in) {
        read_pod(in, block_bytes_);
        if (!valid_block_bytes(block_bytes_)) throw std::invalid_argument("unsupported size of the blocks of the BWT");
        visit([&](auto &ranks) { ranks.read(in); });
    }
};

/// One cache line of bits: a count of 16 bits and 496 bits.
struct alignas(64) bit_block {
    /// The low 16 bits of words[0] hold the number of ones before this block,
    /// relative to the superblock; bit t is stored at bit position t + 16 of the line.
    uint64_t words[8];
};
static_assert(sizeof(bit_block) == 64, "bit_block must be one cache line");

/// Rank over a plain bitvector; 1.03 bits per bit.
class bit_rank {
  public:
    static constexpr size_t HEADER_BITS = 16;
    static constexpr size_t BLOCK_BITS = 512 - HEADER_BITS;
    /// 128 blocks hold 63,488 bits, so the relative counts fit into 16 bits.
    static constexpr size_t SUPER_BLOCKS = 128;

  private:
    static constexpr uint64_t HEADER_MASK = (1ULL << HEADER_BITS) - 1;

    size_t size_ = 0;
    aligned_array<bit_block> blocks_;
    std::vector<uint64_t> supers_;

  public:
    bit_rank() = default;
    explicit bit_rank(size_t size) : size_(size), blocks_(size / BLOCK_BITS + 1) {}

    size_t size() const { return size_; }

    size_t bytes() const { return blocks_.size() * sizeof(bit_block) + supers_.size() * sizeof(uint64_t); }

    /// Set bit i to one; valid only before finalize().
    void set(size_t i) {
        size_t p = i % BLOCK_BITS + HEADER_BITS;
        blocks_[i / BLOCK_BITS].words[p / 64] |= 1ULL << (p % 64);
    }

    /// Compute the rank counts; to be called once after all bits are set.
    void finalize(unsigned threads = 1) {
        // As in quad_rank: the superblocks are counted independently and their totals summed up afterwards.
        size_t count = (blocks_.size() + SUPER_BLOCKS - 1) / SUPER_BLOCKS;
        supers_.assign(blocks_.size() / SUPER_BLOCKS + 1, 0);
        std::vector<uint64_t> sums(count, 0);
        size_t chunk = std::max<size_t>(1, parallel_chunk_items / (SUPER_BLOCKS * BLOCK_BITS));
        parallel_chunks(count, threads, chunk, [&](size_t begin, size_t end, unsigned) {
            for (size_t s = begin; s < end; ++s) {
                uint64_t total = 0;
                size_t last = std::min(blocks_.size(), (s + 1) * SUPER_BLOCKS);
                for (size_t b = s * SUPER_BLOCKS; b < last; ++b) {
                    uint64_t ones = 0;
                    for (int w = 0; w < 8; ++w) ones += __builtin_popcountll(blocks_[b].words[w]);
                    blocks_[b].words[0] |= total;
                    total += ones;
                }
                sums[s] = total;
            }
        });
        uint64_t total = 0;
        for (size_t s = 0; s < count; ++s) {
            supers_[s] = total;
            total += sums[s];
        }
    }

    inline bool operator[](size_t i) const {
        size_t p = i % BLOCK_BITS + HEADER_BITS;
        return (blocks_[i / BLOCK_BITS].words[p / 64] >> (p % 64)) & 1;
    }

    inline void prefetch(size_t i) const { __builtin_prefetch(&blocks_[i / BLOCK_BITS]); }

    /// Number of ones among the first i bits.
    inline size_t rank(size_t i) const {
        size_t b = i / BLOCK_BITS;
        int bits = (int)(i % BLOCK_BITS + HEADER_BITS);
        const bit_block &block = blocks_[b];
        size_t result = supers_[b / SUPER_BLOCKS] + (block.words[0] & HEADER_MASK);
        result += __builtin_popcountll(block.words[0] & ~HEADER_MASK & low_bits(bits));
        for (int w = 1; w < 8; ++w) {
            result += __builtin_popcountll(block.words[w] & low_bits(bits - 64 * w));
        }
        return result;
    }

    /// First step of select(r).
    inline select_state select_prepare(size_t r) const {
        // The last superblock with at most r ones before it; the guess assumes them evenly spread.
        const size_t count = (blocks_.size() + SUPER_BLOCKS - 1) / SUPER_BLOCKS;
        const uint64_t before_last = supers_[count - 1];
        size_t guess = before_last > 0 ? (size_t)((unsigned __int128)r * (count - 1) / before_last) : 0;
        size_t s = first_above(1, count, guess + 1, [&](size_t x) { return supers_[x] > r; }) - 1;
        select_state state;
        state.rank = r - supers_[s];
        state.first = s * SUPER_BLOCKS;
        state.last = std::min(blocks_.size(), state.first + SUPER_BLOCKS);
        size_t blocks = state.last - state.first;
        size_t in_super = s + 1 < count ? supers_[s + 1] - supers_[s] : blocks * BLOCK_BITS / 2;
        state.guess = state.first + std::min(blocks - 1, state.rank * blocks / std::max<size_t>(1, in_super));
        __builtin_prefetch(&blocks_[state.guess]);
        return state;
    }

    /// Second step of select(r).
    inline size_t select_finish(const select_state &state) const {
        // The last block of the superblock with at most that many ones before it; the rank of the state
        // is relative to the superblock.
        size_t b = first_above(state.first + 1, state.last, state.guess + 1,
                               [&](size_t x) { return (blocks_[x].words[0] & HEADER_MASK) > state.rank; }) - 1;
        const bit_block &block = blocks_[b];
        size_t q = state.rank - (block.words[0] & HEADER_MASK);
        for (int w = 0; w < 8; ++w) {
            uint64_t word = w == 0 ? block.words[0] & ~HEADER_MASK : block.words[w];
            size_t ones = __builtin_popcountll(word);
            if (q < ones) return b * BLOCK_BITS + 64 * w + select_in_word(word, q) - HEADER_BITS;
            q -= ones;
        }
        return size_;
    }

    /// Position of the one number r (from 0); there must be more than r ones.
    inline size_t select(size_t r) const { return select_finish(select_prepare(r)); }

    void write(std::ostream &out) const {
        write_pod(out, (uint64_t)size_);
        write_pod(out, (uint64_t)blocks_.size());
        out.write(reinterpret_cast<const char *>(blocks_.data()), blocks_.size() * sizeof(bit_block));
        write_pod(out, (uint64_t)supers_.size());
        out.write(reinterpret_cast<const char *>(supers_.data()), supers_.size() * sizeof(uint64_t));
    }

    void read(std::istream &in) {
        uint64_t size, block_count, super_count;
        read_pod(in, size);
        read_pod(in, block_count);
        size_ = size;
        blocks_.resize(block_count);
        in.read(reinterpret_cast<char *>(blocks_.data()), block_count * sizeof(bit_block));
        read_pod(in, super_count);
        supers_.resize(super_count);
        in.read(reinterpret_cast<char *>(supers_.data()), super_count * sizeof(uint64_t));
    }
};

/// One cache line of a Tunstall-coded bitvector: two counts of 16 bits and 120 codewords of 4 bits.
struct alignas(64) tunstall_block {
    /// Number of bits before this block modulo 2^16, and number of ones before it relative to the superblock.
    uint16_t start;
    uint16_t ones;
    /// Codeword c is the low (c even) or high (c odd) half of codes[c / 2].
    uint8_t codes[60];
};
static_assert(sizeof(tunstall_block) == 64, "tunstall_block must be one cache line");

/// Rank over a bitvector compressed by a variable-to-fixed (Tunstall) code.
///
/// The SA-transformed mask is close to a sequence of independent bits, so only its density can be
/// exploited. The bits are parsed into 16 strings chosen for the density (the leaves of the Tunstall
/// tree), each coded by 4 bits. Codewords of a fixed size make the decoding of a block cheap (the
/// lengths of all its codewords are looked up at once) and waste no space at the end of a block. A
/// block then holds a variable number of bits, so the block of a position is predicted by interpolating
/// in a directory, sparse enough to stay in the cache, and corrected by the counts stored in the blocks.
class tunstall_bits {
  public:
    static constexpr size_t CODE_BYTES = 60;
    static constexpr size_t MAX_LENGTH = 15;
    /// 32 blocks hold at most 57,600 bits, so the relative counts of ones fit into 16 bits.
    static constexpr size_t SUPER_BLOCKS = 32;
    /// The directory has one entry per 16,384 positions (about 28 blocks). The blocks considered for a
    /// position start less than 2^15 bits from it, so the starts of the blocks are needed modulo 2^16 only.
    static constexpr size_t BUCKET_SHIFT = 14;
    /// Fractional bits of the block coordinates in the directory.
    static constexpr size_t FRACTION_SHIFT = 16;

  private:
    /// The strings of the code: their lengths, bits (the first bit is the lowest one) and numbers of ones.
    uint8_t length_[16];
    uint16_t pattern_[16];
    uint8_t ones_[16];
    /// Total lengths and numbers of ones of the two codewords of a byte.
    uint8_t byte_length_[256];
    uint8_t byte_ones_[256];

    size_t size_ = 0;
    aligned_array<tunstall_block> blocks_;
    /// Number of ones before each superblock.
    std::vector<uint64_t> supers_;
    /// The block coordinate of the first position of each bucket, and of the end: the index of the block
    /// with the position plus the fraction of the block before it (FRACTION_SHIFT bits).
    std::vector<uint64_t> directory_;

    void derive_tables() {
        for (int c = 0; c < 16; ++c) ones_[c] = __builtin_popcount(pattern_[c]);
        for (int b = 0; b < 256; ++b) {
            byte_length_[b] = length_[b & 15] + length_[b >> 4];
            byte_ones_[b] = ones_[b & 15] + ones_[b >> 4];
        }
    }

    /// The Tunstall code of a source of independent bits: the most probable string is extended 14 times.
    void make_code(double density) {
        struct leaf {
            uint16_t pattern;
            uint8_t length;
            double probability;
        };
        std::vector<leaf> leaves = {{0, 1, 1 - density}, {1, 1, density}};
        while (leaves.size() < 16) {
            size_t best = leaves.size();
            for (size_t x = 0; x < leaves.size(); ++x) {
                if (leaves[x].length == MAX_LENGTH) continue;
                if (best == leaves.size() || leaves[x].probability > leaves[best].probability) best = x;
            }
            leaf parent = leaves[best];
            leaves[best] = {parent.pattern, (uint8_t)(parent.length + 1), parent.probability * (1 - density)};
            leaves.push_back({(uint16_t)(parent.pattern | (1u << parent.length)), (uint8_t)(parent.length + 1),
                              parent.probability * density});
        }
        for (int c = 0; c < 16; ++c) {
            length_[c] = leaves[c].length;
            pattern_[c] = leaves[c].pattern;
        }
        derive_tables();
    }

    /// The block coordinate of position i interpolated between the ends of its bucket.
    inline uint64_t coordinate(size_t i) const {
        const uint64_t *entry = &directory_[i >> BUCKET_SHIFT];
        return entry[0] + (((entry[1] - entry[0]) * (i & ((1ULL << BUCKET_SHIFT) - 1))) >> BUCKET_SHIFT);
    }

    inline size_t predict(size_t i) const { return coordinate(i) >> FRACTION_SHIFT; }

    /// The bit at the offset of the block, or -1 if the block ends before the offset.
    /// With [counting], also the number of ones of the block before the offset.
    template <bool counting>
    inline int decode(const tunstall_block &block, size_t offset, size_t &ones_before) const {
        uint8_t codes;
#ifdef FMSI2_NEON
        // Lengths of the two codewords of each byte, for the whole block at once; the counts have none.
        alignas(16) static const uint8_t SKIP_COUNTS[16] = {0, 0, 0, 0, 255, 255, 255, 255,
                                                            255, 255, 255, 255, 255, 255, 255, 255};
        const uint8x16_t length_table = vld1q_u8(length_), ones_table = vld1q_u8(ones_);
        alignas(16) uint8_t lengths[64], ones[64];
        uint32x2_t length_sums[4], ones_sums[4];
        for (int v = 0; v < 4; ++v) {
            uint8x16_t x = vld1q_u8(reinterpret_cast<const uint8_t *>(&block) + 16 * v);
            uint8x16_t low = vandq_u8(x, vdupq_n_u8(15)), high = vshrq_n_u8(x, 4);
            uint8x16_t sum = vaddq_u8(vqtbl1q_u8(length_table, low), vqtbl1q_u8(length_table, high));
            if (v == 0) sum = vandq_u8(sum, vld1q_u8(SKIP_COUNTS));
            vst1q_u8(lengths + 16 * v, sum);
            // Totals of the two groups of 8 bytes.
            length_sums[v] = vmovn_u64(vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(sum))));
            if (counting) {
                sum = vaddq_u8(vqtbl1q_u8(ones_table, low), vqtbl1q_u8(ones_table, high));
                if (v == 0) sum = vandq_u8(sum, vld1q_u8(SKIP_COUNTS));
                vst1q_u8(ones + 16 * v, sum);
                ones_sums[v] = vmovn_u64(vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(sum))));
            }
        }
        // The ends of the 8 groups; the group with the offset is the number of groups ending at or before it.
        uint16x8_t ends = vcombine_u16(vmovn_u32(vcombine_u32(length_sums[0], length_sums[1])),
                                       vmovn_u32(vcombine_u32(length_sums[2], length_sums[3])));
        const uint16x8_t zero = vdupq_n_u16(0);
        ends = vaddq_u16(ends, vextq_u16(zero, ends, 7));
        ends = vaddq_u16(ends, vextq_u16(zero, ends, 6));
        ends = vaddq_u16(ends, vextq_u16(zero, ends, 4));
        uint16x8_t passed = vcleq_u16(ends, vdupq_n_u16((uint16_t)std::min(offset, (size_t)0xFFFF)));
        size_t group = vaddvq_u16(vshrq_n_u16(passed, 15));
        if (group == 8) return -1;
        offset -= vmaxvq_u16(vandq_u16(ends, passed));
        // Likewise the byte within the group; a multiplication gives the prefix sums of its 8 lengths.
        uint64_t group_lengths;
        memcpy(&group_lengths, lengths + 8 * group, 8);
        uint8x8_t prefix = vcreate_u8(group_lengths * 0x0101010101010101ULL);
        uint8x8_t skipped = vcle_u8(prefix, vdup_n_u8((uint8_t)offset));
        size_t byte = vaddv_u8(vshr_n_u8(skipped, 7));
        offset -= vmaxv_u8(vand_u8(prefix, skipped));
        codes = reinterpret_cast<const uint8_t *>(&block)[8 * group + byte];
        if (counting) {
            uint16x8_t group_ones = vcombine_u16(vmovn_u32(vcombine_u32(ones_sums[0], ones_sums[1])),
                                                 vmovn_u32(vcombine_u32(ones_sums[2], ones_sums[3])));
            ones_before = vaddvq_u16(vandq_u16(group_ones, passed)) +
                          vaddv_u8(vand_u8(vld1_u8(ones + 8 * group), skipped));
        }
#else
        size_t byte = 0;
        ones_before = 0;
        while (byte < CODE_BYTES && offset >= byte_length_[block.codes[byte]]) {
            offset -= byte_length_[block.codes[byte]];
            ones_before += byte_ones_[block.codes[byte++]];
        }
        if (byte == CODE_BYTES) return -1;
        codes = block.codes[byte];
#endif
        if (offset >= length_[codes & 15]) {
            offset -= length_[codes & 15];
            if (counting) ones_before += ones_[codes & 15];
            codes >>= 4;
        }
        if (counting) ones_before += __builtin_popcount(pattern_[codes & 15] & ((1u << offset) - 1));
        return (pattern_[codes & 15] >> offset) & 1;
    }

    /// The block with position i, searched from block b of its bucket; the bit and, with [counting],
    /// the number of ones of the block before it.
    template <bool counting>
    inline size_t locate(size_t i, size_t b, int &bit, size_t &ones_before) const {
        for (;;) {
            // The distance from the start of the block, known modulo 2^16.
            int16_t offset = (int16_t)((uint16_t)i - blocks_[b].start);
            if (offset < 0) {
                --b;
                continue;
            }
            bit = decode<counting>(blocks_[b], offset, ones_before);
            if (bit >= 0) return b;
            ++b;
        }
    }

  public:
    tunstall_bits() = default;

    size_t size() const { return size_; }

    size_t bytes() const {
        return blocks_.size() * sizeof(tunstall_block) + supers_.size() * sizeof(uint64_t) +
               directory_.size() * sizeof(uint64_t);
    }

    /// Build from the bits get(0), ..., get(size - 1).
    template <typename getter>
    void build(size_t size, getter get) {
        size_ = size;
        size_t total_ones = 0;
        for (size_t i = 0; i < size; ++i) total_ones += get(i);
        make_code(size ? (double)total_ones / size : 0.5);

        // The parsing tree: a child is an inner node, or a leaf given as the negated codeword minus one.
        std::vector<std::array<int, 2>> tree(1, {0, 0});
        for (int c = 0; c < 16; ++c) {
            int node = 0;
            for (int t = 0; t < length_[c]; ++t) {
                int bit = (pattern_[c] >> t) & 1;
                if (t + 1 == length_[c]) {
                    tree[node][bit] = -c - 1;
                } else {
                    if (tree[node][bit] == 0) {
                        tree[node][bit] = tree.size();
                        tree.push_back({0, 0});
                    }
                    node = tree[node][bit];
                }
            }
        }
        // Parse the bits, followed by zeros until position [size] is covered and the last block is full.
        auto parse = [&](auto emit) {
            size_t codes = 0, i = 0;
            int node = 0;
            while (i <= size || node != 0 || codes % (2 * CODE_BYTES) != 0) {
                int next = tree[node][i < size ? (int)get(i) : 0];
                ++i;
                if (next < 0) {
                    emit(codes++, -next - 1);
                    node = 0;
                } else {
                    node = next;
                }
            }
            return codes;
        };
        size_t codes = parse([](size_t, int) {});
        blocks_.resize(codes / (2 * CODE_BYTES));
        supers_.assign(blocks_.size() / SUPER_BLOCKS + 1, 0);
        // The directory: the buckets starting within each block, and the end of the last bucket.
        directory_.assign((size >> BUCKET_SHIFT) + 2, (uint64_t)blocks_.size() << FRACTION_SHIFT);
        size_t bits = 0, ones = 0, first = 0;
        parse([&](size_t code_index, int code) {
            size_t b = code_index / (2 * CODE_BYTES), slot = code_index % (2 * CODE_BYTES);
            if (slot == 0) {
                if (b % SUPER_BLOCKS == 0) supers_[b / SUPER_BLOCKS] = ones;
                blocks_[b].start = (uint16_t)bits;
                blocks_[b].ones = (uint16_t)(ones - supers_[b / SUPER_BLOCKS]);
                first = bits;
            }
            blocks_[b].codes[slot / 2] |= code << (4 * (slot % 2));
            bits += length_[code];
            ones += ones_[code];
            if (slot + 1 == 2 * CODE_BYTES) {
                for (size_t bucket = (first + ((1ULL << BUCKET_SHIFT) - 1)) >> BUCKET_SHIFT;
                     bucket + 1 < directory_.size() && (bucket << BUCKET_SHIFT) < bits; ++bucket) {
                    directory_[bucket] = ((uint64_t)b << FRACTION_SHIFT) +
                                         ((((bucket << BUCKET_SHIFT) - first) << FRACTION_SHIFT) / (bits - first));
                }
            }
        });
    }

    /// Prefetch the directory entries needed by prefetch(i).
    inline void prefetch_directory(size_t i) const { __builtin_prefetch(&directory_[i >> BUCKET_SHIFT]); }

    /// Prefetch the block predicted for position i and return it, to be passed to get().
    inline size_t prefetch(size_t i) const {
        uint64_t position = coordinate(i);
        const tunstall_block *block = &blocks_[position >> FRACTION_SHIFT];
        __builtin_prefetch(block);
        // The prediction is off by one mostly for positions close to an end of the block.
        uint64_t fraction = position & ((1ULL << FRACTION_SHIFT) - 1);
        if (fraction < (1ULL << (FRACTION_SHIFT - 2))) __builtin_prefetch(block - 1);
        if (fraction > (3ULL << (FRACTION_SHIFT - 2))) __builtin_prefetch(block + 1);
        return position >> FRACTION_SHIFT;
    }

    /// Bit i, given a block of its bucket (the one returned by prefetch(i)).
    inline bool get(size_t i, size_t block) const {
        int bit;
        size_t ones_before;
        locate<false>(i, block, bit, ones_before);
        return bit;
    }

    /// Whether the block of position i is predicted right (for tests and statistics).
    bool predicted(size_t i) const {
        size_t b = predict(i), ones_before;
        int16_t offset = (int16_t)((uint16_t)i - blocks_[b].start);
        return offset >= 0 && decode<false>(blocks_[b], offset, ones_before) >= 0;
    }

    inline bool operator[](size_t i) const { return get(i, predict(i)); }

    /// Number of ones among the first i bits.
    inline size_t rank(size_t i) const {
        int bit;
        size_t ones_before;
        size_t b = locate<true>(i, predict(i), bit, ones_before);
        return supers_[b / SUPER_BLOCKS] + blocks_[b].ones + ones_before;
    }

    void write(std::ostream &out) const {
        write_pod(out, (uint64_t)size_);
        out.write(reinterpret_cast<const char *>(length_), sizeof(length_));
        out.write(reinterpret_cast<const char *>(pattern_), sizeof(pattern_));
        write_pod(out, (uint64_t)blocks_.size());
        out.write(reinterpret_cast<const char *>(blocks_.data()), blocks_.size() * sizeof(tunstall_block));
        write_pod(out, (uint64_t)supers_.size());
        out.write(reinterpret_cast<const char *>(supers_.data()), supers_.size() * sizeof(uint64_t));
        write_pod(out, (uint64_t)directory_.size());
        out.write(reinterpret_cast<const char *>(directory_.data()), directory_.size() * sizeof(uint64_t));
    }

    void read(std::istream &in) {
        uint64_t size, count;
        read_pod(in, size);
        size_ = size;
        in.read(reinterpret_cast<char *>(length_), sizeof(length_));
        in.read(reinterpret_cast<char *>(pattern_), sizeof(pattern_));
        derive_tables();
        read_pod(in, count);
        blocks_.resize(count);
        in.read(reinterpret_cast<char *>(blocks_.data()), count * sizeof(tunstall_block));
        read_pod(in, count);
        supers_.resize(count);
        in.read(reinterpret_cast<char *>(supers_.data()), count * sizeof(uint64_t));
        read_pod(in, count);
        directory_.resize(count);
        in.read(reinterpret_cast<char *>(directory_.data()), count * sizeof(uint64_t));
    }
};

} // namespace fmsi2
