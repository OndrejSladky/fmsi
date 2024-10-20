#pragma once

#include <string>
#include <iostream>
#include <cstdint>

#include "models.h"

#ifdef LARGE_KMERS
    typedef __uint128_t kmer_t;
    constexpr int KMER_T_SIZE = 128;
#else
    typedef uint64_t kmer_t;
    constexpr int KMER_T_SIZE = 64;
#endif

/// Convert the given basic nucleotide to int so it can be used for indexing in AC.
/// If non-existing nucleotide is given, return -1.
int NucleotideToInt (char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case 'a': return 0;
        case 'c': return 1;
        case 'g': return 2;
        case 't': return 3;
        default: return -1;
    }
}

/// Convert the given k-mer to its representation as a number.
kmer_t KMerToNumber(const KMer &kMer) {
    kmer_t ret = 0;
    for (char c : kMer.value) {
        ret <<= 2;
        ret |= NucleotideToInt(c);
    }
    return ret;
}

/// Compute the prefix of size d of the given k-mer.
kmer_t BitPrefix(kmer_t kMer, int k, int d) {
    return kMer >> ((k - d) << kmer_t(1));
}

/// Compute the suffix of size d of the given k-mer.
kmer_t BitSuffix(kmer_t kMer, int d) {
    return kMer & ((kmer_t(1) << (d << kmer_t(1))) - kmer_t(1));
}

/// Checkered mask. cmask<uint16_t, 1> is every other bit on
/// (0x55). cmask<uint16_t,2> is two bits one, two bits off (0x33). Etc.
/// Copyright: Jellyfish GPL-3.0
template<typename U, int len, int l = sizeof(U) * 8 / (2 * len)>
struct cmask {
    static const U v =
            (cmask<U, len, l - 1>::v << (2 * len)) | (((U)1 << len) - 1);
};
template<typename U, int len>
struct cmask<U, len, 0> {
    static const U v = 0;
};

/// Compute the reverse complement of a word.
/// Copyright: Jellyfish GPL-3.0
inline kmer_t word_reverse_complement(kmer_t w) {
    typedef kmer_t U;
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
#ifdef LARGE_KMERS
    w = ((w >> 32) & cmask<U, 32>::v) | ((w & cmask<U, 32>::v) << 32);
    w = ( w >> 64                   ) | ( w                    << 64);
#else
    w = ( w >> 32                   ) | ( w                    << 32);
#endif
    return ((U)-1) - w;
}

/// Compute the reverse complement of the given k-mer.
kmer_t ReverseComplement(kmer_t kMer, int k) {
    return (((kmer_t)word_reverse_complement(kMer)) >> (KMER_T_SIZE - (k << kmer_t(1)))) & ((kmer_t(1) << (k << kmer_t (1))) - kmer_t(1));
}

/// Return the complementary nucleotide for the given one.
char ComplementaryNucleotide(const char nucleotide) {
    if (nucleotide == 'A') return 'T';
    else if (nucleotide == 'T') return 'A';
    else if (nucleotide == 'G') return 'C';
    else if (nucleotide == 'C') return 'G';
    throw std::invalid_argument("cannot find complementary nucleotide for letter " + std::string(1, nucleotide));
}

/// Compute the reverse complement of the given k-mer.
KMer ReverseComplement(const KMer &kMer) {
    KMer ans;
    ans.value.reserve(kMer.length());
    for (int i = (int)kMer.length() - 1; i >= 0; --i) {
        ans.value += ComplementaryNucleotide(kMer.value[i]);
    }
    return ans;
}

const char letters[4] {'A', 'C', 'G', 'T'};

/// Get the uppercase version of the character.
inline char UpperToLower(char c) {
    return c - 'A' + 'a';
}

/// Return the index-th nucleotide from the encoded k-mer.
inline char NucleotideAtIndex(kmer_t encoded, int k, int index) {
    return letters[(encoded >> ((k - index - kmer_t(1)) << kmer_t(1))) & kmer_t(3)];
}

/// Convert the encoded KMer representation to string.
std::string NumberToKMer(kmer_t encoded, int length) {
    std::string ret(length, 'N');
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret[length - i -1] = letters[encoded & 3];
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}
