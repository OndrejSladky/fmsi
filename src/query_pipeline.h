#pragma once

// Query pipeline of the FMSI 2 prototype.
//
// The queries flow through three stages, as in the pipeline of minimap2 (kt_pipeline):
//  - a reader parses the query file and cuts it into blocks of whole sequences; a block is
//    closed once it holds a given number of bases, so the work is split by bases and not by
//    the number of sequences;
//  - workers take the blocks one by one as they become free (dynamic scheduling), search
//    their k-mers and format their output;
//  - a writer prints the outputs in the order of the input.
// With one thread, the three stages run in turn on that thread. With more threads, the reader
// and the writer get a thread each, so parsing and printing overlap with the search.
//
// The blocks come from a pool of a fixed size and return to it once printed, which bounds
// the memory and keeps the reader from running ahead of the search.

#include <charconv>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "fms_index2.h"
#include "parser.h"

namespace fmsi2 {

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
    /// Names and sequences of the records, concatenated.
    std::string names, sequences;
    std::vector<record> records;
    std::string output;
};

/// Cuts the query file into blocks.
struct block_reader {
    kseq_t *seq;
    /// A block is closed once it has this many bases; a longer sequence is never split.
    size_t block_bases;
    /// The record held by seq does not belong to any block yet.
    bool pending = false;
    bool too_long = false;

    /// Fill the block with the next sequences; false at the end of the input.
    bool read(query_block &b) {
        b.names.clear();
        b.sequences.clear();
        b.records.clear();
        while (b.sequences.size() < block_bases) {
            if (!pending && kseq_read(seq) < 0) break;
            pending = false;
            // Positions in a block are 32-bit.
            if (b.sequences.size() + seq->seq.l >= UINT32_MAX) {
                if (b.records.empty()) {
                    too_long = true;
                } else {
                    pending = true;
                }
                break;
            }
            b.records.push_back({(uint32_t)b.names.size(), (uint32_t)seq->name.l, (uint32_t)b.sequences.size(),
                                 (uint32_t)seq->seq.l});
            b.names.append(seq->name.s, seq->name.l);
            b.sequences.append(seq->seq.s, seq->seq.l);
        }
        return !b.records.empty();
    }
};

struct query_stats {
    lookup_stats lookups;
    /// Seconds spent in the stages; those of the workers are summed over the threads.
    double read = 0, scan = 0, search = 0, format = 0, write = 0;
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

    void process(const index &ix, query_block &b, bool streaming, unsigned batch_size) {
        const size_t k = ix.k;
        const char *sequences = b.sequences.data();
        auto start = std::chrono::steady_clock::now();

        // Stretches of valid characters.
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
    }
};

/// Query the k-mers of all sequences of the file and print a line per sequence: its name, a tabulator and
/// the presence (0/1) of its k-mers, or their comma-separated hashes (-1 if absent) if [hashes].
/// Returns false if a sequence is too long.
template <bool hashes>
bool run_queries(const index &ix, kseq_t *seq, FILE *out, bool streaming, unsigned batch_size, unsigned threads,
                 size_t block_bases, query_stats &stats) {
    block_reader reader{seq, block_bases};
    std::vector<block_worker<hashes>> workers(threads);

    if (threads == 1) {
        query_block b;
        for (;;) {
            auto start = std::chrono::steady_clock::now();
            bool more = reader.read(b);
            stats.read += seconds_since(start);
            if (!more) break;
            workers[0].process(ix, b, streaming, batch_size);
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
                    workers[t].process(ix, *b, streaming, batch_size);
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
        stats.scan += worker.stats.scan;
        stats.search += worker.stats.search;
        stats.format += worker.stats.format;
    }
    return !reader.too_long;
}

} // namespace fmsi2
