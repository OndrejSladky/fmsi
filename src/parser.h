#pragma once

#include "kseq.h"
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <zlib.h>
#include <vector>
#include "kmers.h"

KSEQ_INIT(gzFile, gzread)

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


/// Obtain k based on the mask convention of k-1 trailing zeros.
int infer_k(std::string ms) {
    int k = 1;
    while(!is_upper(ms[ms.size() - k])) {
        k++;
    }
    return k;
}


/// Load masked superstring from the fasta file in the mask-cased format.
std::string read_masked_superstring(std::string fn) {
    gzFile fp = OpenFile(fn);
    kseq_t *seq = kseq_init(fp);
  int64_t l = kseq_read(seq);
  if (l < 0) {
    throw std::invalid_argument("Error reading the fasta file. The fasta file should contain a single entry - the masked superstring.");
  }
  std::string ret = std::string (seq->seq.s, seq->seq.s + l);
  if ((l = kseq_read(seq)) >= 0) {
    std::cerr << "Warning: The fasta file contains more than one entry. Only the first entry will be used." << std::endl;
  }
  kseq_destroy(seq);
  gzclose(fp);
  return ret;
}

size_t next_invalid_character_or_end(char* sequence, size_t length) {
    for (size_t i = 0; i < length; ++i) {
        if (nucleotideToInt[(uint8_t)sequence[i]] == 4) {
            return i;
        }
    }
    return length;
}
