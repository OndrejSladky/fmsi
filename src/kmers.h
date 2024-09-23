#pragma once

static const uint8_t nucleotideToInt[] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
static const char complementaryNucleotide[] = {
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'A', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 't', 'N', 'g',  'N', 'N', 'N', 'c',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'a', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

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
template <typename U>
inline U word_reverse_complement(U w) {
    w = ((w >> 2)  & cmask<U, 2 >::v) | ((w & cmask<U, 2 >::v) << 2);
    w = ((w >> 4)  & cmask<U, 4 >::v) | ((w & cmask<U, 4 >::v) << 4);
    w = ((w >> 8)  & cmask<U, 8 >::v) | ((w & cmask<U, 8 >::v) << 8);
    w = ((w >> 16) & cmask<U, 16>::v) | ((w & cmask<U, 16>::v) << 16);
    w = ( w >> 32                   ) | ( w                    << 32);
    return ((U)-1) - w;
}

/// Compute the reverse complement of the given k-mer.
template <typename kmer_t>
kmer_t ReverseComplement(kmer_t kMer, int k) {
    return (((kmer_t)word_reverse_complement(kMer)) >> ((sizeof(kMer)<<3) - (k << 1))) & ((kmer_t(1) << (k << 1)) - kmer_t(1));
}

/// Reverse complement a string in place.
void ReverseComplementString(char* s, size_t length) {
    for (size_t i = 0; i < length / 2; ++i) {
        char c = s[i];
        s[i] = complementaryNucleotide[(uint8_t)s[length - i - 1]];
        s[length - i - 1] = complementaryNucleotide[(uint8_t)c];
    }
    if (length % 2) {
        s[length / 2] = complementaryNucleotide[(uint8_t)s[length / 2]];
    }
}


inline bool is_upper(char c) {
    return c >= 'A' && c <= 'Z';
}
