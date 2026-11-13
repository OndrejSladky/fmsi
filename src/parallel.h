#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <vector>

/// Number of items of a loop handed to a thread at a time.
inline size_t parallel_chunk_items = 1 << 20;

/// The size of a chunk of a parallel loop whose chunk boundaries must be multiples of [align]:
/// threads then never write to the same word if the items are bits or symbols packed into words of [align] items.
inline size_t parallel_chunk(size_t align) { return (parallel_chunk_items + align - 1) / align * align; }

/// Call f(begin, end, thread) for consecutive chunks of [0, n) of the given size; the threads take the
/// chunks as they become free.
template <typename function>
void parallel_chunks(size_t n, unsigned threads, size_t chunk, function f) {
    size_t chunks = (n + chunk - 1) / chunk;
    threads = (unsigned)std::max<size_t>(1, std::min<size_t>(threads, chunks));
    std::atomic<size_t> next{0};
    auto work = [&](unsigned thread) {
        for (size_t c; (c = next.fetch_add(1)) < chunks;) f(c * chunk, std::min(n, (c + 1) * chunk), thread);
    };
    std::vector<std::thread> pool;
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto &thread : pool) thread.join();
}
