#pragma once

// Rank structures of the FMSI 2 prototype.
//
// Both structures store the data and the rank counts of a block in the same
// 64-byte cache line, so that a rank query touches a single line (plus a small
// superblock array that stays in the cache). This follows the one-line rank of
// SBWT-squeeze (Diseth et al., 2026) and QuadRank (Groot Koerkamp, 2026).

#include <algorithm>
#include <array>
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
