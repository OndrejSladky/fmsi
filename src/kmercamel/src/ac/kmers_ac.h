#pragma once
#include <string>
#include <vector>
#include "../kmers.h"

struct KMer {
    std::string value;
    size_t length() const{ return value.size(); }
};

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
kmer64_t KMerToNumber(const KMer &kMer) {
    kmer64_t ret = 0;
    for (char c : kMer.value) {
        ret <<= 2;
        ret |= NucleotideToInt(c);
    }
    return ret;
}

