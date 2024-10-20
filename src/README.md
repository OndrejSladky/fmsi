# FMSI code documentation

This file contains a very high-level overview of how the code of FMSI is structured and of algorithms used.

## Code structure

- `main.cpp` contains logic for parsing command line arguments and calling low level functions.
- `fms_index.h` contains all the main index logic, including searching both single and streaming queries, index construction, original masked superstring retrieval and storing and loading the index.
- `parser.h` contains a wrapper around `kseq.h` which parses FASTA files.
- `kmers.h` contains some very basic functions for handling k-mers and strings.
- `function.h` contains definition of some *demasking functions* which are used for experimental support for set operations over k-mers.
- `compact.h` contains a compaction algorithm (which uses KmerCamel🐫's local greedy) used as a subprocedure in the experimental support of set operations.

## Index overview

The implemented index consists of the following parts:

- BWT (fields `ac_gt`, `ac`, and `gt`) requiring 2 bits per character (1 bit for `ac_gt` and 1 bit for both `ac` and `gt`) and associated ranks (0.125 bits per character)
- SA-trasnformed mask (`sa_transformed_mask`) stored with RRR.
- Counts (`counts`), where `counts[c]` is the number of lexicographically smaller characters than `c`.
- Dollar position (`dollar_position`) is the position of the end of the string character in BWT. In BWT this is otherwise stored as `A`.
- The $k$LCP array (`klcp`) stored as plain bit-vector.
- The predictor of strands the queries are from (`predictor`).

The implemented single queries (`query_kmers_single`) and stremaing queries (`query_kmers_streaming`) distinguish three cases:
- `all` which is for masked superstrings that maximize the number of ones (the fasted mode which is supposed to be used).
- `or` is for arbitrary masked superstrings (should be used if it is for some reason undesirable to maximize ones).
- `general` which includes all functions from `function.h` and is meant to be used only within the experimental framework of set operations.


