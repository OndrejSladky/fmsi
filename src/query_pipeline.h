#pragma once

// Query pipeline of the FMSI 2 prototype.
//
// The queries flow through three stages, as in the pipeline of minimap2 (kt_pipeline):
//  - a reader cuts the query file into blocks of whole sequences; a block is closed once it
//    holds a given number of bytes, so the work is split by bases and not by the number of
//    sequences;
//  - workers take the blocks one by one as they become free (dynamic scheduling), search
//    their k-mers and format their output;
//  - a writer prints the outputs in the order of the input.
// With one thread, the three stages run in turn on that thread. With more threads, the reader
// and the writer get a thread each, so reading and printing overlap with the search.
//
// The reader does as little as it can, as it is a single thread. For FASTA it does not parse at
// all: it cuts the bytes of the file before a line starting with '>' and the workers parse their
// blocks, with the same parser (kseq) as for a whole file, so that the sequences are the same.
// This cannot be done within the quality of a FASTQ record, where a line may start with any
// character; the reader therefore parses the sequences itself if the file starts with '@' and
// from the first line starting with '+' on.
//
// The blocks come from a pool of a fixed size and return to it once printed, which bounds
// the memory and keeps the reader from running ahead of the search.

#include <atomic>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "fms_index2.h"
#include "parser.h"

namespace fmsi2 {

/// Bytes for the parser: those in memory first, then the rest of a file, if any.
struct byte_source {
    const char *data = nullptr;
    size_t size = 0, position = 0;
    gzFile file = nullptr;
};

inline int byte_source_read(byte_source *source, void *buffer, unsigned length) {
    if (source->position < source->size) {
        size_t count = std::min<size_t>(length, source->size - source->position);
        memcpy(buffer, source->data + source->position, count);
        source->position += count;
        return (int)count;
    }
    return source->file != nullptr ? gzread(source->file, buffer, length) : 0;
}

/// The parser of FASTA and FASTQ of the rest of FMSI (kseq), reading from a byte_source.
namespace memseq {
KSEQ_INIT2(static inline, byte_source *, byte_source_read)
}

/// A queue of a bounded size passing items between threads.
template <typename T>
class channel {
    std::mutex mutex;
    std::condition_variable not_empty, not_full;
    std::deque<T> items;
    size_t capacity;
    bool closed = false;

  public:
    explicit channel(size_t capacity) : capacity(capacity) {}

    void push(T item) {
        std::unique_lock<std::mutex> lock(mutex);
        not_full.wait(lock, [&] { return items.size() < capacity; });
        items.push_back(std::move(item));
        not_empty.notify_one();
    }

    /// Wait for an item; false once the channel is closed and empty.
    bool pop(T &item) {
        std::unique_lock<std::mutex> lock(mutex);
        not_empty.wait(lock, [&] { return !items.empty() || closed; });
        if (items.empty()) return false;
        item = std::move(items.front());
        items.pop_front();
        not_full.notify_one();
        return true;
    }

    /// No more items will be pushed.
    void close() {
        std::lock_guard<std::mutex> lock(mutex);
        closed = true;
        not_empty.notify_all();
    }
};

/// A block of query sequences: the unit of work passed from the reader through a worker to the writer.
struct query_block {
    struct record {
        uint32_t name, name_length, offset, length;
    };
    /// Position of the block in the input.
    size_t id = 0;
    /// The text of whole FASTA records, if the block is still to be parsed.
    std::string raw;
    bool parsed = false;
    /// Names and sequences of the records, concatenated.
    std::string names, sequences;
    std::vector<record> records;
    std::string output;

    void clear() {
        raw.clear();
        names.clear();
        sequences.clear();
        records.clear();
    }

    /// Append a record; false if the block would get too long for its 32-bit positions.
    bool append(const kstring_t &name, const kstring_t &sequence) {
        if (sequences.size() + sequence.l >= UINT32_MAX) return false;
        records.push_back({(uint32_t)names.size(), (uint32_t)name.l, (uint32_t)sequences.size(), (uint32_t)sequence.l});
        names.append(name.s, name.l);
        sequences.append(sequence.s, sequence.l);
        return true;
    }
};

/// Parse the text of a block; false if a sequence is too long.
///
/// The text has no line starting with '+' (the reader sees to it) and ends at the end of the file or
/// with the end of a line. For such a text, this gives the names and sequences that kseq gives, with
/// one search for the end of each line and one copy of each sequence:
///  - a record starts at the first '>' or '@' of the text and then at each line starting with one;
///  - its name ends at the first white space, the rest of the line is a comment;
///  - its sequence is made of the following lines, without the empty ones and without a carriage
///    return at the end of a line (unless it is all there is of the sequence, or the last byte of the text).
inline bool parse_block(query_block &b) {
    const char *p = b.raw.data(), *end = p + b.raw.size();
    while (p < end && *p != '>' && *p != '@') ++p;
    bool fits = true;
    while (p < end && fits) {
        const char *name = ++p;
        // A header at the very end of the file, with nothing after its first character, is no record.
        if (name == end) break;
        while (p < end && !isspace((unsigned char)*p)) ++p;
        size_t name_length = p - name;
        if (p < end && *p++ != '\n') {
            const char *line_end = (const char *)memchr(p, '\n', end - p);
            p = line_end != nullptr ? line_end + 1 : end;
        }
        size_t offset = b.sequences.size();
        while (p < end && *p != '>' && *p != '@') {
            const char *line = p++;
            if (*line == '\n') continue;
            bool more = p < end;
            const char *line_end = (const char *)memchr(p, '\n', end - p);
            if (line_end == nullptr) line_end = end;
            b.sequences.append(line, line_end - line);
            p = line_end < end ? line_end + 1 : end;
            if (more && b.sequences.size() - offset > 1 && b.sequences.back() == '\r') b.sequences.pop_back();
        }
        fits = b.sequences.size() < UINT32_MAX;
        if (fits) {
            b.records.push_back({(uint32_t)b.names.size(), (uint32_t)name_length, (uint32_t)offset,
                                 (uint32_t)(b.sequences.size() - offset)});
            b.names.append(name, name_length);
        }
    }
    b.parsed = true;
    return fits;
}

/// Cuts the query file into blocks.
class block_reader {
    gzFile file;
    /// A block is closed once it has this many bytes; a longer sequence is never split.
    size_t block_bytes;
    /// Number of bytes read from the file at a time.
    size_t read_step;
    /// Bytes read from the file and not yet handed out; they start at the start of a record.
    std::string buffer;
    /// Bytes of the buffer already searched for a line starting with '+'.
    size_t checked = 0;
    bool at_end = false;
    enum class phase { undecided, cutting, parsing, finished };
    phase step = phase::undecided;
    /// For parsing: the buffer followed by the rest of the file.
    byte_source source;
    memseq::kseq_t *seq = nullptr;
    /// The record held by seq does not belong to any block yet.
    bool pending = false;

    /// Read from the file until the buffer has at least the given size.
    void fill(size_t size) {
        while (!at_end && buffer.size() < size) {
            size_t old = buffer.size();
            buffer.resize(old + read_step);
            int count = gzread(file, &buffer[old], read_step);
            if (count <= 0) {
                count = 0;
                at_end = true;
            }
            buffer.resize(old + count);
        }
    }

    /// Whether a line of the bytes not checked yet starts with '+'.
    bool plus_line() {
        size_t i = checked;
        checked = buffer.size();
        while ((i = buffer.find('+', i)) != std::string::npos) {
            if (i == 0 || buffer[i - 1] == '\n') return true;
            ++i;
        }
        return false;
    }

    void start_parsing() {
        source = byte_source{buffer.data(), buffer.size(), 0, at_end ? nullptr : file};
        seq = memseq::kseq_init(&source);
        step = phase::parsing;
    }

    /// The next whole records of the buffer, as they are.
    bool cut(query_block &b) {
        if (buffer.size() < block_bytes) fill(block_bytes);
        // In the quality of a FASTQ record, a line starting with '>' does not start a record.
        if (plus_line()) {
            start_parsing();
            return parse(b);
        }
        if (buffer.empty()) {
            step = phase::finished;
            return false;
        }
        size_t end = buffer.size();
        if (!at_end || buffer.size() > block_bytes) {
            // The last record starting within the block; else the first one after it, read as far as needed.
            size_t start = block_bytes >= 2 ? buffer.rfind("\n>", block_bytes - 2) : std::string::npos;
            size_t from = block_bytes >= 2 ? block_bytes - 1 : 0;
            while (start == std::string::npos) {
                start = buffer.find("\n>", from);
                if (start != std::string::npos || at_end) break;
                from = buffer.size() - 1;
                fill(buffer.size() + block_bytes);
                if (plus_line()) {
                    start_parsing();
                    return parse(b);
                }
            }
            // Without a further record, the block takes all that has been read.
            end = start != std::string::npos ? start + 1 : buffer.size();
        }
        b.raw.assign(buffer, 0, end);
        b.parsed = false;
        buffer.erase(0, end);
        checked -= end;
        return true;
    }

    /// The next records, parsed.
    bool parse(query_block &b) {
        b.parsed = true;
        while (b.sequences.size() < block_bytes && b.records.size() < block_bytes) {
            // As when a whole file is parsed, the first error ends the input.
            if (!pending && memseq::kseq_read(seq) < 0) {
                step = phase::finished;
                break;
            }
            pending = !b.append(seq->name, seq->seq);
            if (pending) {
                if (b.records.empty()) {
                    too_long = true;
                    step = phase::finished;
                }
                break;
            }
        }
        return !b.records.empty();
    }

  public:
    std::atomic<bool> too_long{false};

    block_reader(gzFile file, size_t block_bytes, size_t read_step = 1 << 20)
        : file(file), block_bytes(std::max<size_t>(block_bytes, 1)), read_step(std::max<size_t>(read_step, 1)) {}
    block_reader(const block_reader &) = delete;
    ~block_reader() { memseq::kseq_destroy(seq); }

    /// Fill the block with the next sequences; false at the end of the input.
    bool read(query_block &b) {
        b.clear();
        if (step == phase::undecided) {
            // The parser skips everything before the first '>' or '@'.
            fill(block_bytes);
            size_t first = buffer.find_first_of(">@");
            while (first == std::string::npos && !at_end) {
                buffer.clear();
                fill(block_bytes);
                first = buffer.find_first_of(">@");
            }
            if (first == std::string::npos) {
                step = phase::finished;
            } else if (buffer[first] == '@') {
                start_parsing();
            } else {
                step = phase::cutting;
            }
        }
        switch (step) {
        case phase::cutting:
            return cut(b);
        case phase::parsing:
            return parse(b);
        default:
            return false;
        }
    }
};

struct query_stats {
    lookup_stats lookups;
    /// Seconds spent in the stages; those of the workers are summed over the threads.
    double read = 0, parse = 0, scan = 0, search = 0, format = 0, write = 0;
};

inline double seconds_since(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

/// Searches the blocks and formats their output; one per thread.
template <bool hashes>
struct block_worker {
    typedef typename std::conditional<hashes, int64_t, int8_t>::type out_t;
    std::vector<uint32_t> offsets;
    std::vector<segment> segments;
    std::vector<out_t> results;
    query_stats stats;

    /// False if a sequence of the block is too long.
    bool process(const index &ix, query_block &b, bool streaming, unsigned batch_size) {
        const size_t k = ix.k;
        auto start = std::chrono::steady_clock::now();
        bool fits = b.parsed || parse_block(b);
        stats.parse += seconds_since(start);
        const char *sequences = b.sequences.data();

        // Stretches of valid characters.
        start = std::chrono::steady_clock::now();
        offsets.clear();
        segments.clear();
        for (const auto &r : b.records) {
            size_t valid = 0;
            for (size_t i = 0; i <= r.length; ++i) {
                if (i < r.length && nucleotideToInt[(uint8_t)sequences[r.offset + i]] != 4) {
                    ++valid;
                    if (!streaming && valid >= k) offsets.push_back(r.offset + i + 1 - k);
                } else {
                    if (streaming && valid >= k) segments.push_back({(uint32_t)(r.offset + i - valid), (uint32_t)valid});
                    valid = 0;
                }
            }
        }
        // K-mers with invalid characters are reported as absent.
        results.assign(b.sequences.size(), hashes ? -1 : 0);
        stats.scan += seconds_since(start);

        start = std::chrono::steady_clock::now();
        if (streaming) {
            lookup_streaming<hashes, out_t>(ix, sequences, segments.data(), segments.size(), results.data(), batch_size,
                                            stats.lookups);
        } else {
            lookup_batch<hashes, out_t>(ix, sequences, offsets.data(), offsets.size(), results.data(), batch_size,
                                        stats.lookups);
        }
        stats.search += seconds_since(start);

        start = std::chrono::steady_clock::now();
        std::string &out = b.output;
        out.clear();
        for (const auto &r : b.records) {
            out.append(b.names, r.name, r.name_length);
            out.push_back('\t');
            size_t count = r.length >= k ? r.length - k + 1 : 0;
            if (hashes) {
                char number[24];
                for (size_t i = 0; i < count; ++i) {
                    if (i > 0) out.push_back(',');
                    char *end = std::to_chars(number, number + sizeof(number), (int64_t)results[r.offset + i]).ptr;
                    out.append(number, end - number);
                }
            } else {
                size_t at = out.size();
                out.resize(at + count);
                for (size_t i = 0; i < count; ++i) out[at + i] = (char)('0' + results[r.offset + i]);
            }
            out.push_back('\n');
        }
        stats.format += seconds_since(start);
        return fits;
    }
};

/// Query the k-mers of all sequences of the file and print a line per sequence: its name, a tabulator and
/// the presence (0/1) of its k-mers, or their comma-separated hashes (-1 if absent) if [hashes].
/// Returns false if a sequence is too long.
template <bool hashes>
bool run_queries(const index &ix, gzFile file, FILE *out, bool streaming, unsigned batch_size, unsigned threads,
                 size_t block_bases, query_stats &stats) {
    block_reader reader(file, block_bases);
    std::vector<block_worker<hashes>> workers(threads);

    if (threads == 1) {
        query_block b;
        for (;;) {
            auto start = std::chrono::steady_clock::now();
            bool more = reader.read(b);
            stats.read += seconds_since(start);
            if (!more) break;
            if (!workers[0].process(ix, b, streaming, batch_size)) reader.too_long = true;
            start = std::chrono::steady_clock::now();
            fwrite(b.output.data(), 1, b.output.size(), out);
            stats.write += seconds_since(start);
        }
    } else {
        // Blocks in flight; the reader waits for a printed block once all are in use.
        const size_t pool_size = 4 * (size_t)threads;
        std::vector<query_block> pool(pool_size);
        channel<query_block *> unused(pool_size), to_search(pool_size), to_print(pool_size);
        for (auto &b : pool) unused.push(&b);

        std::thread reading([&] {
            query_block *b;
            for (size_t id = 0; unused.pop(b); ++id) {
                auto start = std::chrono::steady_clock::now();
                bool more = reader.read(*b);
                stats.read += seconds_since(start);
                if (!more) break;
                b->id = id;
                to_search.push(b);
            }
            to_search.close();
        });
        std::vector<std::thread> searching;
        for (unsigned t = 0; t < threads; ++t) {
            searching.emplace_back([&, t] {
                query_block *b;
                while (to_search.pop(b)) {
                    if (!workers[t].process(ix, *b, streaming, batch_size)) reader.too_long = true;
                    to_print.push(b);
                }
            });
        }
        std::thread writing([&] {
            // Blocks finished ahead of their turn.
            std::map<size_t, query_block *> waiting;
            size_t next = 0;
            query_block *b;
            while (to_print.pop(b)) {
                waiting.emplace(b->id, b);
                while (!waiting.empty() && waiting.begin()->first == next) {
                    b = waiting.begin()->second;
                    waiting.erase(waiting.begin());
                    ++next;
                    auto start = std::chrono::steady_clock::now();
                    fwrite(b->output.data(), 1, b->output.size(), out);
                    stats.write += seconds_since(start);
                    unused.push(b);
                }
            }
        });
        reading.join();
        for (auto &thread : searching) thread.join();
        to_print.close();
        writing.join();
    }

    for (const auto &worker : workers) {
        stats.lookups.add(worker.stats.lookups);
        stats.parse += worker.stats.parse;
        stats.scan += worker.stats.scan;
        stats.search += worker.stats.search;
        stats.format += worker.stats.format;
    }
    return !reader.too_long;
}

} // namespace fmsi2
