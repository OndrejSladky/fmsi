#pragma once

#include <vector>
#include <list>

#include "kmers.h"
#include "khash.h"


#define kh_int128_hash_func(key) kh_int64_hash_func((khint64_t)((key)>>65^(key)^(key)<<21))
#define kh_int128_hash_equal(a, b) ((a) == (b))
#define kh_int256_hash_func(key) kh_int128_hash_func((__uint128_t)((key)>>129^(key)^(key)<<35))
#define kh_int256_hash_equal(a, b) ((a) == (b))

#define KHASH_MAP_INIT_INT128(name, khval_t)								\
	KHASH_INIT(name, __uint128_t, khval_t, 1, kh_int128_hash_func, kh_int128_hash_equal)

#define KHASH_SET_INIT_INT128(name)										\
    KHASH_INIT(name, __uint128_t, char, 0, kh_int128_hash_func, kh_int128_hash_equal)

#define KHASH_MAP_INIT_INT256(name, khval_t)								\
	KHASH_INIT(name, uint256_t, khval_t, 1, kh_int256_hash_func, kh_int256_hash_equal)

#define KHASH_SET_INIT_INT256(name)										\
	KHASH_INIT(name, uint256_t, char, 0, kh_int256_hash_func, kh_int256_hash_equal)

// Use 128-bit integers for extra large k-mers to allow for larger k.
KHASH_SET_INIT_INT256(S256)
KHASH_MAP_INIT_INT256(P256, size_t)
// Use 128-bit integers for large k-mers to allow for larger k.
KHASH_SET_INIT_INT128(S128)
KHASH_MAP_INIT_INT128(P128, size_t)
// Use 64-bits integers for small k-mers for faster operations and less memory usage.
KHASH_SET_INIT_INT64(S64)
KHASH_MAP_INIT_INT64(P64, size_t)

#define INIT_KHASH_WRAPPER(type) \
    struct kmer_dict##type##_t { \
        inline kh_S##type##_t *kh_init_set() { \
            return kh_init_S##type(); \
        }                         \
        inline khint_t kh_get_from_set(kh_S##type##_t *set, kmer##type##_t key) { \
            return kh_get_S##type(set, key); \
        }                        \
        inline khint_t kh_put_to_set(kh_S##type##_t *set, kmer##type##_t key, int *ret) { \
            return kh_put_S##type(set, key, ret); \
        }                        \
        inline void kh_del_from_set(kh_S##type##_t *set, khint_t key) { \
            kh_del_S##type(set, key); \
        }                        \
        inline void kh_destroy_set(kh_S##type##_t *set) { \
            kh_destroy_S##type(set); \
        }                        \
        inline kh_P##type##_t *kh_init_map() { \
            return kh_init_P##type(); \
        }                         \
        inline khint_t kh_get_from_map(kh_P##type##_t *map, kmer##type##_t key) { \
            return kh_get_P##type(map, key); \
        }                        \
        inline khint_t kh_put_to_map(kh_P##type##_t *map, kmer##type##_t key, int *ret) { \
            return kh_put_P##type(map, key, ret); \
        }                        \
        inline void kh_del_from_map(kh_P##type##_t *map, khint_t key) { \
            kh_del_P##type(map, key); \
        }                        \
        inline void kh_destroy_map(kh_P##type##_t *map) { \
            kh_destroy_P##type(map); \
        }                        \
        inline void kh_clear_map(kh_P##type##_t *map) { \
            kh_clear_P##type(map); \
        }                        \
        inline void kh_resize_map(kh_P##type##_t *map, khint_t size) { \
            kh_resize_P##type(map, size); \
        }                        \
    };

INIT_KHASH_WRAPPER(64)
INIT_KHASH_WRAPPER(128)
INIT_KHASH_WRAPPER(256)

/// Determine whether the k-mer or its reverse complement is present.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
bool containsKMer(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t kMer, int k, bool complements) {
    bool ret = wrapper.kh_get_from_set(kMers, kMer) != kh_end(kMers);
    if (complements) ret |= wrapper.kh_get_from_set(kMers, ReverseComplement(kMer, k )) != kh_end(kMers);
    return ret;
}

/// Remove the k-mer and its reverse complement.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void eraseKMer(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t kMer, int k, bool complements) {
    auto key = wrapper.kh_get_from_set(kMers, kMer);
    if (key != kh_end(kMers)) {
        wrapper.kh_del_from_set(kMers, key);
    }
    if (complements) {
        kmer_t reverseComplement = ReverseComplement(kMer, k);
        key = wrapper.kh_get_from_set(kMers, reverseComplement);
        if (key != kh_end(kMers)) wrapper.kh_del_from_set(kMers, key);
    }
}

/// Return the next k-mer in the k-mer set and update the index.
template <typename kmer_t, typename kh_S_t>
kmer_t nextKMer(kh_S_t *kMers, [[maybe_unused]] kmer_t _, size_t &lastIndex) {
    for (size_t i = kh_begin(kMers) + lastIndex; i != kh_end(kMers); ++i, ++lastIndex) {
        if (!kh_exist(kMers, i)) continue;
        return kh_key(kMers, i);
    }
    // No more k-mers.
    return -1;
}

/// Construct a vector of the k-mer set in an arbitrary order.
template <typename kmer_t, typename kh_S_t>
std::vector<kmer_t> kMersToVec(kh_S_t *kMers, [[maybe_unused]] kmer_t _) {
    std::vector<kmer_t> res(kh_size(kMers));
    size_t index = 0;
    for (auto i = kh_begin(kMers); i != kh_end(kMers); ++i) {
        if (!kh_exist(kMers, i)) continue;
        res[index++] = kh_key(kMers, i);
    }
    return res;
}

/// Add an interval with given index to the given k-mer.
///
/// [intervalsForKMer] store the intervals and [intervals] maps the k-mer to the index in [intervalsForKMer].
template <typename kmer_t, typename kh_P_t, typename kh_wrapper_t>
bool appendInterval(kh_P_t *intervals, kh_wrapper_t wrapper, std::vector<std::list<size_t>> &intervalsForKMer, kmer_t kMer, size_t index, int k, bool complements) {
    if (complements) kMer = std::min(kMer, ReverseComplement(kMer, k));
    auto key = wrapper.kh_get_from_map(intervals, kMer);
    if (key == kh_end(intervals)) {
        int ret;
        wrapper.kh_put_to_map(intervals, kMer, &ret);
        key = wrapper.kh_get_from_map(intervals, kMer);
        kh_value(intervals, key) = intervalsForKMer.size();
        intervalsForKMer.emplace_back(std::list<size_t>());
    }
    key = wrapper.kh_get_from_map(intervals, kMer);
    auto position = kh_value(intervals, key);
    if (intervalsForKMer[position].empty() || intervalsForKMer[position].back() != index) {
        intervalsForKMer[position].push_back(index);
        return true;
    }
    return false;
}




