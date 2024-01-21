#pragma once

#include "kseq.h"
#include "mask.h"
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

typedef struct {
  mask_t mask;
  std::string superstring;
} masked_superstring_t;

inline bool _isupper(const char c) { return c >= 'A' && c <= 'Z'; }

inline char _toupper(const char c) { return _isupper(c) ? c : (c - 'a' + 'A'); }

/// Load masked superstring from the fasta file in the mask-cased format.
masked_superstring_t read_masked_superstring(std::string fn) {
  gzFile fp;
  kseq_t *seq;
  masked_superstring_t ret;

  fp = gzopen(fn.c_str(), "r");
  if (fp == 0) {
    throw std::invalid_argument("couldn't open file " + fn);
  }

  seq = kseq_init(fp);

  int l;
  while ((l = kseq_read(seq)) >= 0) {
    for (int i = 0; i < l; ++i) {
      char c = seq->seq.s[i];
      ret.mask.push_back(_isupper(c));
      ret.superstring.push_back(_toupper(c));
    }
  }
  kseq_destroy(seq);
  gzclose(fp);
  return ret;
}

/// Return the complement of the given nucleotide.
inline char _complement(char c) {
  switch (c) {
  case 'A':
    return 'T';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'T':
    return 'A';
  default:
    return 'N';
  }
}

/// Return the reverse complement of the given sequence.
std::string reverse_complement(std::string s) {
  std::string res;
  for (size_t i = 0; i < s.size(); ++i) {
    res.push_back(_complement(s[s.size() - 1 - i]));
  }
  return res;
}

/// Write the given superstring to file.
void write_superstring(std::string fn, std::string superstring) {
  std::ofstream of;
  of.open(fn, std::ios::out);
  of.write(superstring.c_str(), superstring.size());
  of.close();
}

/// Get the path where to store the mask.
inline std::string compute_mask_path(std::string &fn, int k, bool include_k) {
  return fn + (include_k ? ".k" + std::to_string(k) : "") + ".mask";
}
