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


/// Load masked superstring from the fasta file in the mask-cased format.
std::string read_masked_superstring(std::string fn) {
    gzFile fp = OpenFile(fn);
    kseq_t *seq = kseq_init(fp);
    std::vector<char> ret;
  int64_t l;
  while ((l = kseq_read(seq)) >= 0) {
    for (int64_t i = 0; i < l; ++i) {
      ret.push_back(seq->seq.s[i]);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);
  return std::string(ret.begin(), ret.end());
}

size_t next_invalid_character_or_end(char* sequence, size_t length) {
    for (size_t i = 0; i < length; ++i) {
        if (nucleotideToInt[(uint8_t)sequence[i]] == 4) {
            return i;
        }
    }
    return length;
}
