#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "kmers.h"
#include "khash_utils.h"


/// Fill the k-mer dictionary with k-mers from the given sequence.
/// If complements is true, always add the canonical k-mers.
/// If case_sensitive is true, add the k-mer only if it starts with an upper case letter.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void AddKMers(kh_S_t *kMers, kh_wrapper_t wrapper, [[maybe_unused]] kmer_t _, size_t sequence_length,
              const char* sequence, int64_t k, bool complements, bool case_sensitive = false) {
    int64_t currentLength = 0;
    kmer_t currentKMer = 0, reverseComplement = 0;
    kmer_t cases = 0;
    kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
    kmer_t shift = 2 * (k - 1);
    for (size_t i = 0; i < sequence_length; ++i) {
        auto data = nucleotideToInt[(uint8_t)sequence[i]];
        if (data >= 4) {
            // Restart if "N"-like nucleotide.
            currentKMer = reverseComplement = 0;
            currentLength = 0;
            continue;
        }
        currentKMer = ((currentKMer << 2) | data) & mask;
        reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
        // K-mer is present if it is upper case or case-insensitive.
        cases = (cases | (!case_sensitive || sequence[i] <= 'Z')) << 1;
        if ((++currentLength >= k) && (cases & (kmer_t(1) << k))) {
            // Add the canonical k-mer to the dictionary.
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            int ret;
            wrapper.kh_put_to_set(kMers, canonical, &ret);
        }
    }
}

/// Return a file/stdin for reading.
gzFile OpenFile(std::string &path) {
    FILE *in_stream;
    if(path=="-"){
        in_stream = stdin;
    }
    else {
        in_stream = fopen(path.c_str(), "r");
        if (in_stream == nullptr) {
            throw std::invalid_argument("couldn't open file " + path);
        }
    }
    gzFile fp = gzdopen(fileno(in_stream), "r");
    return fp;
}


/// Load a dictionary of k-mers from a fasta file.
/// If complements is true, add the canonical k-mers.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void ReadKMers(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t _, std::string &path, int k, bool complements,
               bool case_sensitive = false) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        AddKMers(kMers, wrapper, _, seq->seq.l, seq->seq.s, k, complements, case_sensitive);
    }

    kseq_destroy(seq);
    gzclose(fp);
}

/// Read the masked superstring from the given path and return it wrapped as a kseq_t.
kseq_t* ReadMaskedSuperstring(std::string &path) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);
    kseq_read(seq);
    return seq;
}

/// Ensure that the file is at the end.
void AssertEOF(kseq_t *seq, std::string message) {
    if (kseq_read(seq) >= 0) {
        throw std::invalid_argument(message);
    }
}

/// Print the fasta file header.
void WriteName(const int k, std::ostream &of) {
    of << ">superstring ";
    of << "k=" << k << std::endl;
}
