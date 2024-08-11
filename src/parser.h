#pragma once

#include "kseq.h"
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <zlib.h>
#include <vector>

KSEQ_INIT(gzFile, gzread)

/// Load masked superstring from the fasta file in the mask-cased format.
std::string read_masked_superstring(std::string fn) {
  gzFile fp;
  kseq_t *seq;

  fp = gzopen(fn.c_str(), "r");
  if (fp == 0) {
    throw std::invalid_argument("couldn't open file " + fn);
  }

  seq = kseq_init(fp);
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
