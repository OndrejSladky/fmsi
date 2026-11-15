#pragma once

// The mask stored in the order of the text (index files written with fmsi index --compact-file).
//
// In the order of the suffix array the mask is close to a sequence of independent bits, so it compresses only to
// its entropy (0.2 to 1 bit per character). Along the superstring it is a few long runs: the k-mers of a simplitig
// are set and the k-mers that span the junction to the next one are not. Stored as the lengths of its runs it takes
// a fraction of that (0.02 to 0.35 bits per character on the datasets of the paper).
//
// The index needs the mask in the order of the suffix array, so it is permuted when the file is loaded, with the
// BWT alone: the LF-mapping goes from the row of the suffix at position p to that of the suffix at p - 1, so a
// walk from the end of the text visits the rows of all positions in turn, and each row gets the mask bit of its
// position. One walk cannot be split among threads or overlapped with others, so the file also stores the rows of
// the suffixes at every SAMPLE_STEP-th position: the walks between two samples are independent, run on all
// threads, and several of them advance in lockstep, each step prefetching the block of the BWT of its next row.
// Every walk must arrive at the row of the next sample, which checks the file.
//
// Loading takes one LF step per character (about 10 ns per character on one thread, where the mask in the order
// of the suffix array takes 2 ns): the price of the smaller file.

#include <atomic>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "parallel.h"

namespace fmsi2 {

/// The distance of the sampled positions, and the number of walks advanced in lockstep by a thread.
static const uint64_t SAMPLE_STEP = 1 << 16;
static const unsigned LOCKSTEP_WALKS = 32;

/// The mask of a masked superstring as its runs: the lengths of the runs of unset and of set characters in turn,
/// starting with a run of unset ones (of length 0 if the first character is set), each as an integer of variable
/// length (7 bits per byte, the low ones first, the high bit of a byte set if another byte follows).
inline std::vector<uint8_t> mask_runs(const std::string &ms) {
    std::vector<uint8_t> runs;
    auto put = [&](uint64_t value) {
        while (value >= 128) {
            runs.push_back((uint8_t)(value | 128));
            value >>= 7;
        }
        runs.push_back((uint8_t)value);
    };
    bool set = false;
    uint64_t length = 0;
    for (char c : ms) {
        const bool upper = c >= 'A' && c <= 'Z';
        if (upper != set) {
            put(length);
            set = upper;
            length = 0;
        }
        ++length;
    }
    put(length);
    return runs;
}

/// The mask of the superstring read from its runs: the bits of any range of positions, found from cursors at
/// every [step]-th position.
class run_mask {
    struct cursor {
        /// The run with the position: its first byte, the number of characters before it, whether it is set.
        uint64_t offset, before;
        bool set;
    };
    const uint8_t *runs_ = nullptr;
    size_t bytes_ = 0;
    uint64_t size_ = 0, step_ = 1;
    std::vector<cursor> cursors_;

    inline uint64_t read(uint64_t &offset) const {
        uint64_t value = 0;
        for (unsigned shift = 0; offset < bytes_; shift += 7) {
            const uint8_t byte = runs_[offset++];
            value |= (uint64_t)(byte & 127) << shift;
            if (byte < 128) return value;
            if (shift > 56) break;
        }
        throw std::invalid_argument("the runs of the mask are damaged");
    }

  public:
    /// The runs must describe exactly [size] characters.
    run_mask(const uint8_t *runs, size_t bytes, uint64_t size, uint64_t step) : runs_(runs), bytes_(bytes), size_(size), step_(step) {
        cursors_.resize(size / step + 1);
        uint64_t offset = 0, at = 0;
        bool set = false;
        while (offset < bytes) {
            const uint64_t start = offset, length = read(offset);
            if (length > size - at) throw std::invalid_argument("the runs of the mask are longer than the superstring");
            // The cursors of the sampled positions within the run.
            for (uint64_t j = (at + step - 1) / step; j * step < at + length; ++j) cursors_[j] = {start, at, set};
            at += length;
            set = !set;
        }
        if (at != size) throw std::invalid_argument("the runs of the mask do not cover the superstring");
        // A position at the very end (size a multiple of the step) has no run; it is never asked for.
    }

    uint64_t size() const { return size_; }

    /// Set the bits of the positions [begin, end) in [out]: position p sets bit p - begin + shift (out is not
    /// cleared).
    void fill(uint64_t begin, uint64_t end, uint64_t *out, uint64_t shift = 0) const {
        if (begin >= end) return;
        const cursor &c = cursors_[begin / step_];
        uint64_t offset = c.offset, at = c.before;
        bool set = c.set;
        while (at < end) {
            const uint64_t length = read(offset);
            if (set) {
                const uint64_t from = std::max(at, begin), to = std::min(at + length, end);
                for (uint64_t p = from; p < to;) {
                    // Whole words where possible.
                    const uint64_t bit = p - begin + shift, room = 64 - (bit & 63), count = std::min(room, to - p);
                    out[bit >> 6] |= (count == 64 ? ~0ULL : ((1ULL << count) - 1)) << (bit & 63);
                    p += count;
                }
            }
            at += length;
            set = !set;
        }
    }
};

/// Bits set from several threads and read in their order, in chunks that are given back once they are read, so
/// that a structure built from the bits can replace them piece by piece.
class chunked_bits {
    static constexpr unsigned CHUNK_SHIFT = 21;
    static constexpr size_t CHUNK_WORDS = (size_t)1 << (CHUNK_SHIFT - 6);
    std::vector<std::unique_ptr<uint64_t[]>> chunks_;
    size_t size_ = 0;

  public:
    explicit chunked_bits(size_t size) : size_(size) {
        chunks_.resize((size >> CHUNK_SHIFT) + 1);
        for (auto &chunk : chunks_) chunk.reset(new uint64_t[CHUNK_WORDS]());
    }

    size_t size() const { return size_; }

    inline void set_atomic(size_t i) {
        __atomic_fetch_or(&chunks_[i >> CHUNK_SHIFT][(i >> 6) & (CHUNK_WORDS - 1)], 1ULL << (i & 63), __ATOMIC_RELAXED);
    }

    inline void prefetch(size_t i) const { __builtin_prefetch(&chunks_[i >> CHUNK_SHIFT][(i >> 6) & (CHUNK_WORDS - 1)], 1); }

    /// The bits in their order, 64 at a time: sink(word, count). With [release], each chunk is freed once read.
    template <typename sink_t>
    void pass(sink_t sink, bool release) {
        for (size_t first = 0; first < size_; first += 64) {
            const size_t chunk = first >> CHUNK_SHIFT;
            sink(chunks_[chunk][(first >> 6) & (CHUNK_WORDS - 1)], (unsigned)std::min<size_t>(64, size_ - first));
            if (release && (first + 64 >= size_ || ((first + 64) >> CHUNK_SHIFT) != chunk)) chunks_[chunk].reset();
        }
    }
};

/// The rows of the suffixes at the positions step, 2 step, ... of the text, as the constructions collect them.
struct text_samples {
    uint64_t step = 0;
    std::vector<uint64_t> rows;
};

/// Permute the mask of the superstring into the order of the suffix array: set(row) is called, from several
/// threads, for every row whose position is set, and prefetch(row) for the rows it may be called for next.
///
/// [ranks] is the rank structure of the BWT, [counts] the cumulative character counts, [dollar] and [dollar2] the
/// rows of the sentinels (dollar2 == NO for one strand). The text is the superstring X of [mask] (strands == 1),
/// or X, a separator and the reverse complement Y of X (strands == 2), in which the k-mer of Y at position
/// n + 1 + q is set iff that of X at n - k - q is. Throws if a walk does not arrive at its sample.
template <typename ranks_t, typename set_t, typename prefetch_t>
void permute_text_mask(const ranks_t &ranks, const uint64_t counts[4], uint64_t dollar, uint64_t dollar2, unsigned strands, unsigned k,
                       const run_mask &mask, const text_samples &samples, unsigned threads, set_t set, prefetch_t prefetch) {
    const uint64_t NO = ~(uint64_t)0;
    const uint64_t n = mask.size(), length = strands == 2 ? 2 * n + 1 : n, step = samples.step;
    if (step == 0 || samples.rows.size() != length / step) throw std::invalid_argument("the samples of the index do not fit its text");
    const uint64_t segments = (length + step - 1) / step;
    const size_t words = (size_t)((step + 63) / 64);

    // The mask of the text positions [a, b): bit p - a.
    auto fill = [&](uint64_t a, uint64_t b, uint64_t *bits, std::vector<uint64_t> &scratch) {
        std::fill(bits, bits + words, 0);
        if (a < n) mask.fill(a, std::min(b, n), bits);
        if (strands == 2 && b > n + 1 && n >= k) {
            // Positions p of Y: the bit of X at i = 2n + 1 - k - p, for i in [0, n - k].
            const uint64_t first = std::max(a, n + 1), top = 2 * n + 1 - k;
            if (first > top) return;
            const uint64_t low = top >= b - 1 ? top - (b - 1) : 0, high = std::min(top - first, n - k);
            if (low > high) return;
            scratch.assign((size_t)((high - low) / 64 + 1), 0);
            mask.fill(low, high + 1, scratch.data());
            for (size_t w = 0; w < scratch.size(); ++w) {
                for (uint64_t word = scratch[w]; word != 0; word &= word - 1) {
                    const uint64_t i = low + 64 * w + (uint64_t)__builtin_ctzll(word), p = top - i;
                    bits[(p - a) >> 6] |= 1ULL << ((p - a) & 63);
                }
            }
        }
    };

    std::atomic<bool> failed{false};
    parallel_chunks(segments, threads, 4 * LOCKSTEP_WALKS, [&](size_t begin, size_t end, unsigned) {
        struct walk {
            uint64_t row, position, start;
            /// A row to be set at the next step of the walk, once its word has been prefetched (NO: none).
            uint64_t pending;
            bool active;
        };
        std::vector<walk> walks(LOCKSTEP_WALKS);
        std::vector<uint64_t> bits(words * LOCKSTEP_WALKS), scratch;
        size_t next = begin, active = 0;
        auto start = [&](walk &w, size_t slot) {
            if (next >= end) {
                w.active = false;
                return;
            }
            const uint64_t a = (uint64_t)next * step, b = std::min(a + step, length);
            ++next;
            fill(a, b, bits.data() + slot * words, scratch);
            // The walk starts at the suffix of position b: the empty one (row 0), or a sampled one.
            w = {b == length ? 0 : samples.rows[b / step - 1], b, a, NO, true};
            ranks.prefetch(w.row);
            ++active;
        };
        for (size_t slot = 0; slot < walks.size(); ++slot) start(walks[slot], slot);
        while (active > 0) {
            for (size_t slot = 0; slot < walks.size(); ++slot) {
                walk &w = walks[slot];
                if (!w.active) continue;
                if (w.pending != NO) {
                    set(w.pending);
                    w.pending = NO;
                }
                // From the row of the suffix at w.position to that of the suffix one position before.
                if (w.row == dollar2) {
                    // The suffix Y: before it is the separator, whose suffix is the second smallest one.
                    w.row = 1;
                    --w.position;
                } else {
                    const uint8_t c = ranks[w.row];
                    const uint64_t sentinels = c == 0 ? (uint64_t)(w.row > dollar) + (uint64_t)(dollar2 != NO && w.row > dollar2) : 0;
                    w.row = counts[c] + ranks.rank_inlined(w.row, c) - sentinels;
                    --w.position;
                    const uint64_t t = w.position - w.start;
                    if ((bits[slot * words + (t >> 6)] >> (t & 63)) & 1) {
                        w.pending = w.row;
                        prefetch(w.row);
                    }
                }
                if (w.position == w.start) {
                    const uint64_t expected = w.start == 0 ? dollar : samples.rows[w.start / step - 1];
                    if (w.row != expected) failed = true;
                    if (w.pending != NO) set(w.pending);
                    --active;
                    start(w, slot);
                } else {
                    ranks.prefetch(w.row);
                }
            }
        }
    });
    if (failed) throw std::invalid_argument("the mask of the index does not fit its BWT: the file is damaged");
}

} // namespace fmsi2
