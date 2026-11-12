#pragma once

// Rank structures of the FMSI 2 prototype.
//
// Both structures store the data and the rank counts of a block in the same
// 64-byte cache line, so that a rank query touches a single line (plus a small
// superblock array that stays in the cache). This follows the one-line rank of
// SBWT-squeeze (Diseth et al., 2026) and QuadRank (Groot Koerkamp, 2026).

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <vector>

// FMSI2_SCALAR selects the portable code also on ARM (for testing).
#if defined(__aarch64__) && !defined(FMSI2_SCALAR)
#include <arm_neon.h>
#define FMSI2_NEON 1
#endif

namespace fmsi2 {

/// A zero-initialized heap array aligned to cache lines.
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
        if (posix_memalign(&memory, 64, size * sizeof(T)) != 0) throw std::bad_alloc();
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

/// One cache line of 2-bit symbols: 4 counts of 16 bits and 224 symbols.
struct alignas(64) quad_block {
    /// Occurrences of each symbol before this block, relative to the superblock.
    uint16_t counts[4];
    /// Symbol t is stored at bits [2t % 64, 2t % 64 + 2) of words[t / 32].
    uint64_t words[7];
};
static_assert(sizeof(quad_block) == 64, "quad_block must be one cache line");

/// Rank over a sequence of 2-bit symbols; 2.29 bits per symbol.
class quad_rank {
  public:
    static constexpr size_t BLOCK_SYMBOLS = 224;
    /// 256 blocks hold 57,344 symbols, so the relative counts fit into 16 bits.
    static constexpr size_t SUPER_BLOCKS = 256;

  private:
    static constexpr uint64_t LOW_BITS_OF_PAIRS = 0x5555555555555555ULL;

    struct super_block {
        uint64_t counts[4];
    };

    /// Byte masks selecting the first r symbols of a block, and no byte of its counts.
    struct limit_table {
        alignas(64) uint8_t masks[BLOCK_SYMBOLS][64];
        limit_table() {
            memset(masks, 0, sizeof(masks));
            for (size_t r = 0; r < BLOCK_SYMBOLS; ++r) {
                for (size_t bit = 0; bit < 2 * r; ++bit) masks[r][8 + bit / 8] |= 1 << (bit % 8);
            }
        }
    };
    static inline const limit_table LIMITS{};

    size_t size_ = 0;
    aligned_array<quad_block> blocks_;
    std::vector<super_block> supers_;

  public:
    quad_rank() = default;
    explicit quad_rank(size_t size) : size_(size), blocks_(size / BLOCK_SYMBOLS + 1) {}

    size_t size() const { return size_; }

    size_t bytes() const {
        return blocks_.size() * sizeof(quad_block) + supers_.size() * sizeof(super_block);
    }

    /// Set symbol i; valid only before finalize() and only once per position.
    void set(size_t i, uint8_t c) {
        size_t r = i % BLOCK_SYMBOLS;
        blocks_[i / BLOCK_SYMBOLS].words[r / 32] |= (uint64_t)c << (2 * (r % 32));
    }

    /// Compute the rank counts; to be called once after all symbols are set.
    void finalize() {
        supers_.assign(blocks_.size() / SUPER_BLOCKS + 1, super_block{{0, 0, 0, 0}});
        uint64_t totals[4] = {0, 0, 0, 0};
        for (size_t b = 0; b < blocks_.size(); ++b) {
            auto &super = supers_[b / SUPER_BLOCKS];
            if (b % SUPER_BLOCKS == 0) {
                for (int c = 0; c < 4; ++c) super.counts[c] = totals[c];
            }
            for (int c = 0; c < 4; ++c) {
                blocks_[b].counts[c] = (uint16_t)(totals[c] - super.counts[c]);
            }
            for (int w = 0; w < 7; ++w) {
                for (int c = 0; c < 4; ++c) {
                    totals[c] += __builtin_popcountll(matches(blocks_[b].words[w], c));
                }
            }
        }
    }

    /// The bits 2t are set for the symbols t of the word that are equal to c.
    static inline uint64_t matches(uint64_t word, uint8_t c) {
        uint64_t x = word ^ (LOW_BITS_OF_PAIRS * (uint64_t)(c ^ 3));
        return x & (x >> 1) & LOW_BITS_OF_PAIRS;
    }

    inline uint8_t operator[](size_t i) const {
        size_t r = i % BLOCK_SYMBOLS;
        return (blocks_[i / BLOCK_SYMBOLS].words[r / 32] >> (2 * (r % 32))) & 3;
    }

    inline void prefetch(size_t i) const { __builtin_prefetch(&blocks_[i / BLOCK_SYMBOLS]); }

    /// Number of occurrences of c among the first i symbols.
    inline size_t rank(size_t i, uint8_t c) const {
        size_t b = i / BLOCK_SYMBOLS;
        size_t r = i % BLOCK_SYMBOLS;
        const quad_block &block = blocks_[b];
        size_t result = supers_[b / SUPER_BLOCKS].counts[c] + block.counts[c];
        // Branch-free on purpose: the position within the block is unpredictable.
#ifdef FMSI2_NEON
        uint8x16_t sum = vdupq_n_u8(0);
        for (int v = 0; v < 4; ++v) {
            sum = vaddq_u8(sum, vcntq_u8(vandq_u8(matches(block, v, c), vld1q_u8(LIMITS.masks[r] + 16 * v))));
        }
        return result + vaddlvq_u8(sum);
#else
        for (int w = 0; w < 7; ++w) {
            result += __builtin_popcountll(matches(block.words[w], c) & low_bits(2 * (int)r - 64 * w));
        }
        return result;
#endif
    }

    /// Ranks of c at the positions i <= j; cheaper than two rank queries if both fall into the same block.
    inline void rank_pair(size_t i, size_t j, uint8_t c, size_t &rank_i, size_t &rank_j) const {
        size_t b = i / BLOCK_SYMBOLS;
        if (b != j / BLOCK_SYMBOLS) {
            rank_i = rank(i, c);
            rank_j = rank(j, c);
            return;
        }
        if (j == i + 1) {
            rank_i = rank(i, c);
            rank_j = rank_i + ((*this)[i] == c);
            return;
        }
        size_t ri = i % BLOCK_SYMBOLS, rj = j % BLOCK_SYMBOLS;
        const quad_block &block = blocks_[b];
        size_t before = supers_[b / SUPER_BLOCKS].counts[c] + block.counts[c];
#ifdef FMSI2_NEON
        uint8x16_t sum_i = vdupq_n_u8(0), sum_j = vdupq_n_u8(0);
        for (int v = 0; v < 4; ++v) {
            uint8x16_t m = matches(block, v, c);
            sum_i = vaddq_u8(sum_i, vcntq_u8(vandq_u8(m, vld1q_u8(LIMITS.masks[ri] + 16 * v))));
            sum_j = vaddq_u8(sum_j, vcntq_u8(vandq_u8(m, vld1q_u8(LIMITS.masks[rj] + 16 * v))));
        }
        rank_i = before + vaddlvq_u8(sum_i);
        rank_j = before + vaddlvq_u8(sum_j);
#else
        rank_i = rank_j = before;
        for (int w = 0; w < 7; ++w) {
            uint64_t m = matches(block.words[w], c);
            rank_i += __builtin_popcountll(m & low_bits(2 * (int)ri - 64 * w));
            rank_j += __builtin_popcountll(m & low_bits(2 * (int)rj - 64 * w));
        }
#endif
    }

#ifdef FMSI2_NEON
    /// The bits 2t of the bytes are set for the symbols equal to c in the v-th quarter of the block.
    static inline uint8x16_t matches(const quad_block &block, int v, uint8_t c) {
        uint8x16_t x = veorq_u8(vld1q_u8(reinterpret_cast<const uint8_t *>(&block) + 16 * v),
                                vdupq_n_u8(0x55 * (c ^ 3)));
        return vandq_u8(vandq_u8(x, vshrq_n_u8(x, 1)), vdupq_n_u8(0x55));
    }
#endif

    void write(std::ostream &out) const {
        write_pod(out, (uint64_t)size_);
        write_pod(out, (uint64_t)blocks_.size());
        out.write(reinterpret_cast<const char *>(blocks_.data()), blocks_.size() * sizeof(quad_block));
        write_pod(out, (uint64_t)supers_.size());
        out.write(reinterpret_cast<const char *>(supers_.data()), supers_.size() * sizeof(super_block));
    }

    void read(std::istream &in) {
        uint64_t size, block_count, super_count;
        read_pod(in, size);
        read_pod(in, block_count);
        size_ = size;
        blocks_.resize(block_count);
        in.read(reinterpret_cast<char *>(blocks_.data()), block_count * sizeof(quad_block));
        read_pod(in, super_count);
        supers_.resize(super_count);
        in.read(reinterpret_cast<char *>(supers_.data()), super_count * sizeof(super_block));
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
    void finalize() {
        supers_.assign(blocks_.size() / SUPER_BLOCKS + 1, 0);
        uint64_t total = 0;
        for (size_t b = 0; b < blocks_.size(); ++b) {
            if (b % SUPER_BLOCKS == 0) supers_[b / SUPER_BLOCKS] = total;
            uint64_t ones = 0;
            for (int w = 0; w < 8; ++w) ones += __builtin_popcountll(blocks_[b].words[w]);
            blocks_[b].words[0] |= total - supers_[b / SUPER_BLOCKS];
            total += ones;
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

} // namespace fmsi2
