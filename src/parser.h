#pragma once

#include "stdexcept"
#include <string>

// TODO remove this
#ifdef USE_MALLOC_WRAPPERS
#undef USE_MALLOC_WRAPPERS
#endif

#include "kseq.h"
#include "mask.h"
#include <stdio.h>
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

typedef struct {
  mask_t mask;
  std::string superstring;
} masked_superstring_t;

inline bool _isupper(const char c) { return c >= 'A' && c <= 'Z'; }

inline char _toupper(const char c) { return _isupper(c) ? c : (c - 'a' + 'A'); }

masked_superstring_t read_masked_superstring(const char *fn) {
  gzFile fp;
  kseq_t *seq;
  masked_superstring_t ret;

  fp = gzopen(fn, "r");
  if (fp == 0) {
    throw std::invalid_argument("couldn't open file");
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

void append_reverse_complement(masked_superstring_t &ms, int k) {
  size_t length = ms.superstring.size();
  for (size_t i = 0; i < length; ++i) {
    ms.superstring.push_back(_complement(ms.superstring[length - 1 - i]));
    ms.mask.push_back(ms.mask[length - k - i]);
  }
}
