#pragma once

// The C++ API of FMSI 2: the supported way to use the index from another program.
//
//   #include "fmsi_api.h"
//
//   fmsi::Index index = fmsi::Index::load("kmers.fmsi");            // an index written by `fmsi index`
//   bool present = index.contains("ACGTACGTACGTACGTACGTACGTACGTACG");
//   std::vector<uint8_t> hits = index.contains_all(read);           // one value per k-mer of the read
//   int64_t hash = index.lookup(kmer);                              // -1 if the k-mer is not in the set
//   std::string kmer = index.access(hash);                          // the inverse of lookup
//
// Two kinds of queries, as in FMSI:
//   - membership (contains...): whether a k-mer is in the set; works with any mask;
//   - dictionary (lookup..., access...): the hash of a k-mer and the k-mer of a hash. With a min-one mask
//     (minimal_hashes()) the hashes are a minimal perfect hash of the set, 0 to number_of_hashes() - 1.
// A k-mer and its reverse complement are the same k-mer of the set. In an index of both strands (preset fast)
// they have different hashes (those of the k-mers as they are), so the hashes range over twice the set.
//
// An Index is immutable once loaded or built: its query methods are const, keep no state between calls and may be
// called from several threads at the same time. The methods for many k-mers (contains_batch, contains_all,
// lookup_batch, lookup_all, access_batch) search several k-mers in lockstep, each step prefetching the memory of
// the next one; they are several times faster per k-mer than the methods for one k-mer, which wait for memory.
//
// Everything else in the sources of FMSI (the namespace fmsi2 and the code of FMSI 0.4) is internal and may change.
// The translation unit that includes this header must be compiled as the Makefile of FMSI compiles main.cpp
// (C++17, -O3 -march=native) and linked with sdsl, zlib and libsais.

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "doubled.h"
#include "fms_index2.h"

namespace fmsi {

/// What is indexed: the masked superstring alone (small), or with its reverse complement (fast: twice the size,
/// about half the query time, as a k-mer is searched as it is only).
enum class Preset { small, fast };

/// The layout of the index in memory, chosen when it is loaded: the one its file recommends (that of its preset),
/// the most compact one (small), or the fastest one (fast: about a sixth larger than small for one strand).
enum class Layout { recommended, small, fast };

struct LoadOptions {
    Layout layout = Layout::recommended;
    /// Load the kLCP array if its file (<index file>.klcp) is there; it speeds up contains_all and lookup_all
    /// on sequences that do not follow the superstring.
    bool klcp = true;
    /// Prepare for access(): the plain mask and the hints of select (1/16 bit per character more). Without it
    /// access works too, more slowly.
    bool access = false;
    /// Threads that compute the rank counts while loading.
    unsigned threads = 1;
    /// Length of the patterns whose first search steps are tabulated; -1 for the default of the layout.
    int table_length = -1;
};

struct BuildOptions {
    Preset preset = Preset::small;
    Layout layout = Layout::recommended;
    /// Also construct the kLCP array.
    bool klcp = false;
    /// As LoadOptions::access.
    bool access = false;
    unsigned threads = 1;
    /// "min-one", "max-one", or "auto" to find the type of the mask out.
    std::string mask_type = "auto";
    int table_length = -1;
};

class Index {
    std::unique_ptr<fmsi2::index> ix_;

    /// Number of k-mers or sequences searched in lockstep.
    static constexpr unsigned BATCH = 32;

    explicit Index(std::unique_ptr<fmsi2::index> ix) : ix_(std::move(ix)) {}

    static fmsi2::layout layout_of(Layout layout, const std::string &preset) {
        if (layout == Layout::recommended) return fmsi2::preset_layout(preset);
        return fmsi2::preset_layout(layout == Layout::fast ? "fast" : "small");
    }

    void prepare(bool access, int table_length) {
        fmsi2::build_prefix_table(*ix_, table_length < 0 ? fmsi2::default_prefix_length(*ix_) : (unsigned)table_length);
        if (access) fmsi2::build_access_support(*ix_);
    }

    static bool valid(const char *kmer, size_t k) {
        for (size_t i = 0; i < k; ++i) {
            if (nucleotideToInt[(uint8_t)kmer[i]] == 4) return false;
        }
        return true;
    }

    /// The k-mers of a sequence: results[i] for the k-mer at i (absent if it has a character other than ACGT).
    template <bool hashes, typename out_t>
    std::vector<out_t> all(std::string_view sequence, bool streaming) const {
        const size_t k = ix_->k, length = sequence.size();
        const size_t count = length >= k ? length - k + 1 : 0;
        if (length >= UINT32_MAX) throw std::invalid_argument("sequences of 4 Gbp and more are not supported");
        std::vector<out_t> results(length, hashes ? -1 : 0);
        std::vector<uint32_t> offsets;
        std::vector<fmsi2::segment> segments;
        size_t valid_run = 0;
        for (size_t i = 0; i <= length; ++i) {
            if (i < length && nucleotideToInt[(uint8_t)sequence[i]] != 4) {
                ++valid_run;
                if (!streaming && valid_run >= k) offsets.push_back((uint32_t)(i + 1 - k));
            } else {
                if (streaming && valid_run >= k) segments.push_back({(uint32_t)(i - valid_run), (uint32_t)valid_run});
                valid_run = 0;
            }
        }
        fmsi2::lookup_stats stats;
        if (streaming) {
            fmsi2::lookup_streaming<hashes, out_t>(*ix_, sequence.data(), segments.data(), segments.size(), results.data(), BATCH, stats);
        } else {
            fmsi2::lookup_batch<hashes, out_t>(*ix_, sequence.data(), offsets.data(), offsets.size(), results.data(), BATCH, stats);
        }
        results.resize(count);
        return results;
    }

    /// [count] k-mers stored one after another: out[i] for the k-mer at kmers + i * k.
    template <bool hashes, typename out_t>
    void batch(const char *kmers, size_t count, out_t *out, unsigned lockstep) const {
        const size_t k = ix_->k;
        // In pieces whose offsets fit into 32 bits; the searches write their results at the offsets of the k-mers.
        const size_t piece = std::max<size_t>(1, std::min<size_t>(1 << 16, (size_t)(UINT32_MAX / 2) / k));
        std::vector<uint32_t> offsets;
        std::vector<out_t> results;
        fmsi2::lookup_stats stats;
        for (size_t first = 0; first < count; first += piece) {
            const size_t n = std::min(piece, count - first);
            const char *sequence = kmers + first * k;
            offsets.clear();
            for (size_t i = 0; i < n; ++i) {
                if (valid(sequence + i * k, k)) offsets.push_back((uint32_t)(i * k));
            }
            results.assign(n * k, hashes ? -1 : 0);
            fmsi2::lookup_batch<hashes, out_t>(*ix_, sequence, offsets.data(), offsets.size(), results.data(), lockstep, stats);
            for (size_t i = 0; i < n; ++i) out[first + i] = results[i * k];
        }
    }

  public:
    Index(Index &&) = default;
    Index &operator=(Index &&) = default;

    /// Load an index file written by `fmsi index` or save() (<name>.fmsi; the name without the extension is
    /// accepted too). Throws std::invalid_argument if the file is missing, damaged, or needs a newer FMSI.
    static Index load(const std::string &path, const LoadOptions &options = LoadOptions()) {
        const std::string file = fmsi2::find_index_file(path);
        if (file.empty()) throw std::invalid_argument("no index file " + path);
        auto ix = std::make_unique<fmsi2::index>();
        fmsi2::layout layout;
        if (options.layout != Layout::recommended) layout = layout_of(options.layout, "");
        if (options.access) layout.mask = (int)fmsi2::mask_kind::plain;
        fmsi2::load(*ix, file, options.threads, layout);
        const std::string klcp = file + fmsi2::KLCP_EXTENSION;
        if (options.klcp && fmsi2::file::file_version(klcp) != 0) fmsi2::load_klcp(*ix, klcp, options.threads);
        Index index(std::move(ix));
        index.prepare(options.access, options.table_length);
        return index;
    }

    /// Build an index in memory from a masked superstring: the superstring with the characters of the k-mers
    /// that are not in the set in lower case, as KmerCamel writes it (without the FASTA header). Construction
    /// takes 5 to 11 bytes of memory per character of the superstring.
    static Index build(std::string masked_superstring, int k, const BuildOptions &options = BuildOptions()) {
        if (k < 1 || (size_t)k > masked_superstring.size()) throw std::invalid_argument("k must be between 1 and the length of the superstring");
        const bool fast = options.preset == Preset::fast;
        const bool with_same_kmer = mbwt_with_same_kmer;
        mbwt_with_same_kmer = options.mask_type == "auto";
        mbwt m = fast ? doubled::compute_mbwt(std::move(masked_superstring), k, options.klcp, options.threads)
                      : k <= 32 ? compute_mbwt<uint64_t>(masked_superstring, k, options.klcp, false, options.threads)
                                : compute_mbwt<__uint128_t>(masked_superstring, k, options.klcp, false, options.threads);
        mbwt_with_same_kmer = with_same_kmer;
        const fmsi2::mask_statistics statistics = fmsi2::decide_mask_type(fmsi2::classify_mask(m), options.mask_type, "");
        fmsi2::layout layout = layout_of(options.layout, fast ? "fast" : "small");
        if (options.access) layout.mask = (int)fmsi2::mask_kind::plain;
        auto ix = std::make_unique<fmsi2::index>();
        fmsi2::build(*ix, m, k, (fmsi2::mask_kind)layout.mask, layout.bwt_block_bytes, options.threads);
        ix->mask_type_name = statistics.type == "unknown" ? "" : statistics.type;
        ix->mask_type_source = statistics.source;
        Index index(std::move(ix));
        index.prepare(options.access, options.table_length);
        return index;
    }

    /// Write the index to a file (.fmsi is added to the name if missing) and its kLCP array, if it has one, to
    /// <file>.klcp. The file recommends the layout of the preset of the index. Returns the name of the file.
    std::string save(const std::string &path, unsigned threads = 1) const {
        const std::string extension = fmsi2::INDEX_EXTENSION;
        const bool has_extension = path.size() >= extension.size() && path.compare(path.size() - extension.size(), extension.size(), extension) == 0;
        const std::string file = has_extension ? path : path + extension;
        const fmsi2::layout recommended = fmsi2::preset_layout(ix_->preset);
        fmsi2::save(*ix_, file, threads, &recommended);
        if (ix_->has_klcp()) fmsi2::save_klcp(*ix_, file + fmsi2::KLCP_EXTENSION, threads);
        return file;
    }

    /// The length of the k-mers.
    unsigned k() const { return ix_->k; }

    /// Whether the reverse complement of the superstring is indexed too (preset fast).
    bool both_strands() const { return ix_->strands == 2; }

    /// Number of hashes: lookup returns values below it, and access takes them. With a min-one mask this is the
    /// number of k-mers of the set; for an index of both strands, the number of k-mers counted in each
    /// orientation (twice the set, less its k-mers that are their own reverse complements).
    uint64_t number_of_hashes() const { return ix_->mask_ones; }

    /// The type of the mask as the index states it: "min-one", "max-one", "min-one,max-one", "other", or "" if unknown.
    const std::string &mask_type() const { return ix_->mask_type_name; }

    /// Whether lookup is a minimal perfect hash of the set: the mask is known to be min-one.
    bool minimal_hashes() const { return ix_->mask_type_name.rfind("min-one", 0) == 0; }

    bool has_klcp() const { return ix_->has_klcp(); }

    /// Bytes of the index in memory.
    size_t bytes() const { return ix_->bytes(); }

    // Membership queries.

    /// Whether the k-mer (k characters) is in the set. One k-mer at a time: for many, use contains_batch.
    bool contains(std::string_view kmer) const {
        if (kmer.size() != ix_->k) throw std::invalid_argument("the k-mer must have k characters");
        uint8_t result = 0;
        batch<false, int8_t>(kmer.data(), 1, reinterpret_cast<int8_t *>(&result), 1);
        return result == 1;
    }

    /// Whether each of [count] k-mers, stored one after another (count * k characters), is in the set: out[i] is
    /// 1 or 0. A k-mer with a character other than ACGT is not in the set.
    void contains_batch(const char *kmers, size_t count, uint8_t *out) const {
        batch<false, int8_t>(kmers, count, reinterpret_cast<int8_t *>(out), BATCH);
    }

    /// Whether each k-mer of the sequence is in the set: one value (1 or 0) per k-mer, the first one for the k-mer
    /// at the start of the sequence. With [streaming], the search of a k-mer continues that of its neighbour (for
    /// reads and longer sequences); without it, the k-mers are searched independently, in lockstep.
    std::vector<uint8_t> contains_all(std::string_view sequence, bool streaming = true) const {
        std::vector<int8_t> found = all<false, int8_t>(sequence, streaming);
        return std::vector<uint8_t>(found.begin(), found.end());
    }

    // Dictionary queries.

    /// The hash of the k-mer, or -1 if it is not in the set.
    int64_t lookup(std::string_view kmer) const {
        if (kmer.size() != ix_->k) throw std::invalid_argument("the k-mer must have k characters");
        int64_t result = -1;
        batch<true, int64_t>(kmer.data(), 1, &result, 1);
        return result;
    }

    /// The hashes of [count] k-mers stored one after another: out[i], -1 for a k-mer that is not in the set.
    void lookup_batch(const char *kmers, size_t count, int64_t *out) const { batch<true, int64_t>(kmers, count, out, BATCH); }

    /// The hash of each k-mer of the sequence (-1 if it is not in the set); see contains_all.
    std::vector<int64_t> lookup_all(std::string_view sequence, bool streaming = true) const {
        return all<true, int64_t>(sequence, streaming);
    }

    /// The k-mer of a hash (below number_of_hashes()), as it occurs in the superstring: the looked-up k-mer or its
    /// reverse complement. The inverse of lookup for a min-one mask.
    std::string access(uint64_t hash) const {
        if (hash >= ix_->mask_ones) throw std::out_of_range("the hash is not below number_of_hashes()");
        std::string kmer(ix_->k, 'A');
        fmsi2::access(*ix_, hash, kmer.data());
        return kmer;
    }

    /// The k-mers of [count] hashes, written one after another to out (count * k characters).
    void access_batch(const uint64_t *hashes, size_t count, char *out) const {
        for (size_t i = 0; i < count; ++i) {
            if (hashes[i] >= ix_->mask_ones) throw std::out_of_range("a hash is not below number_of_hashes()");
        }
        if (ix_->mask_type == fmsi2::mask_kind::plain) {
            fmsi2::access_batch(*ix_, hashes, count, out, BATCH);
        } else {
            // The lockstep access selects on the plain mask; other encodings go one by one.
            for (size_t i = 0; i < count; ++i) fmsi2::access(*ix_, hashes[i], out + i * ix_->k);
        }
    }
};

} // namespace fmsi
