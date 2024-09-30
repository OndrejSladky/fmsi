# FMSI

[![FMSI test](https://github.com/OndrejSladky/fmsi/actions/workflows/ci.yml/badge.svg)](https://github.com/OndrejSladky/fmsi/actions/)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Prerequisities](#prerequisities)
* [Getting started](#getting-started)
* [How to use](#how-to-use)
* [How it works](#how-it-works)
* [How to test](#how-to-test)
* [Issues](#issues)
* [Changelog](#changelog)
* [Licence](#licence)
* [Contact](#contact)

<!-- vim-markdown-toc -->

## Introduction

FMSI is a highly memory-efficient tool (typically 3-5 memory bits / indexed k-mer) for performing membership queries on single $k$-mer sets.
FMSI uses [masked superstrings](https://doi.org/10.1101/2024.03.06.583483) for storing the $k$-mer sets to ensure high compressibility for a wide range of different $k$-mer sets,
and implements FMS-index, a simplification of the FM-index. It supports both streaming and single queries.

It is based on the following papers:

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Function-Assigned Masked Superstrings as a Versatile and Compact Data Type for *k*-Mer Sets.
> *bioRxiv* 2024.03.06.583483, 2024. [https://doi.org/10.1101/2024.03.06.583483](https://doi.org/10.1101/2024.03.06.583483)

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

To construct an index (the `fmsi index` subcommand), FMSI accepts as input (see the `-p` parameter) a masked superstring of the $k$-mer set.
The masked superstring can be computed by [KmerCamel🐫](https://github.com/OndrejSladky/kmercamel).
It then stores the index in files with the same prefix and the `.fmsi` extension.
To query the index (the `fmsi query` subcommand), FMSI accepts a FASTA file with $k$-mers to query (see the `-q` parameter).


Additionally, FMSI experimentally supports set operations using the $f$-masked superstrings framework.
However, this feature is currently only experimental and requires rather significant time and memory.

## Prerequisities

- GCC 4.8+ or equivalent
- Zlib

## Getting started

To compute masked superstrings, download and compile KmerCamel🐫:

```
git clone --recursive https://github.com/OndrejSladky/kmercamel
cd kmercamel && make
```

Compute and optimize a masked superstring:
```
./kmercamel -p ./index.fa -k 31 -c -o ms-not-optimized.fa
./kmercamel optimize -p ms-not-optimized.fa -k 31 -a ones -c -o ms.fa
```

Download and compile FMSI:

```
git clone --recursive https://github.com/OndrejSladky/fmsi
cd prophasm2 && make -j
```


Create an index from a masked superstring:

```
./fmsi index -p ms.fa
```

Query the index (note that the `-O` optimization flag requires optimized masks for the number of ones):

```
./fmsi query -p ms.fa -q query.fa -O
```

## How to use

Construction of the index (get the full list of parameters by running `fmsi index -h`):
```
./fmsi index -p ms.fa        # Index the ms.fa file.
./fmsi index -p ms.fa -s     # Do not construct the kLCP array to save memory on queries if streaming queries are not needed.
./fmsi index -p ms.fa -k 31  # Explicitly provide k to override the automatically inferred value from the mask.
```

Querying the index (get the full list of parameters by running `fmsi query -h`):
```
./fmsi query -p ms.fa -q query.fa -O   # Query the index using optimizations for mask maximizing the number of ones (can lead to incorrect results if the mask is not optimized).
./fmsi query -p ms.fa -q query.fa      # Do not turn on the optimizations.
./fmsi query -p ms.fa                  # Read the queries from stdin.
./fmsi query -p ms.fa -q query.fa -F   # Print results per each entry of the fasta file.
```

### Experimental implementation of set operations

#### Basic usage

If you wish to perform set operations or answer membership queries without the need to understand the
details of the $f$-masked superstring framework, FMSI can manage the details for you.
This can be done in a few simple steps:
1. **Compute a masked superstring.**
   - This can be done by [KmerCamel🐫](https://github.com/OndrejSladky/kmercamel); simply run `kmercamel -c -k 31 -p fasta.fa -o ms.fa` (with appropriate values for `-k` and `-p`).
   - If you obtained the masked superstring in a different way, make sure it minimizes the number of ones; if you're unsure, you can use `kmercamel optimize -c -a zeros -k 31 -p ms_more_ones.fa -o ms.fa`. No need to optimize superstrings directly computed by KmerCamel🐫.
2. **Index the masked superstring** with `./fmsi index -p ms.fa`.
3. **Perform the set operations** you wish with `./fmsi [operation] -k 31 -p ms1.fa -p ms2.fa -r ms.fa`. Possible operations (use the names instead of `[operation]`) are:
   - `union` to compute the union.
   - `inter` to compute the intersection
   - `diff` to compute the set difference
   - `symdiff` to compute the symmetric difference
4. **Query the index** with `./fmsi query -p ms.fa -q query.fa -k 31`.
5. To **get back the underlying masked superstring**, use `./fmsi export -p ms.fa`.

If you use FMSI this way, it ensures that operations in any order and queries to any index are computed correctly,
while keeping the memory usage for queries as low as possible. Furthermore, exported $f$-masked superstrings are always,
or-masked superstrings, which are the default masked superstrings.

The only downside to this approach is that each set operation uses compaction, which is the most time- and memory- consuming
part of the process, which in some use cases might cause slowdowns which are not necessary. If this is your case,
you probably want to stick to the advanced usage, managing the functions and building block methods yourself.

#### Advanced usage

If you wish to manage the operations yourself, the workflow is quite similar to the basic usage, with the following changes:
- Underlying $f$-MS concatenation can be done with `./fmsi merge`. Details on which $f$ should be used is described in Chapter 4 of the [paper](https://doi.org/10.1101/2024.03.06.583483).
- Compaction and change back to the or-masked superstring can be done with `./fmsi compact`.
- If you query the index with different function than or, use the `-f` argument. The same applies to compaction.

### Commands overview

To run the tool, run `./fmsi [command]`

The recognized commands are:

- `index` Creates a BWT based index of the given masked superstring.
- `query` Queries a $k$-mer against an index.
- `union` Performs union on two indexes. Expects or-MS.
- `inter` Performs intersection on two indexes. Expects or-MS.
- `diff`  Performs set difference on two indexes. Expects or-MS.
- `symdiff` Performs symmetric difference on two indexes. Expects or-MS.
- `merge` Merge several indexes.
- `compact` Compacts an index.
- `export` Export the underlying masked superstring.
- `clean` Cleans the files stored for index.
- `-v`    Prints the version of the program.
- `-h`    Prints the help.

Each command has its own set of arguments, which can be displayed by running `./fmsi [command] -h`.

## How it works?

FMSI builds on the representation of *k*-mer sets via [masked superstrings](https://doi.org/10.1101/2024.03.06.583483),
which ensures a generally compact representation of *k*-mers, requiring typically 1-1.4 characters per *k*-mer for genomes and pan-genomes.

FMSI then builds a simplified FM-index, which omits the sampled suffix array, which is the key to the memory efficiency.
It instead computes the SA-transformed mask, which makes it possible to determine the presence of a *k*-mer
without the need to compute the original coordinates.

Additionally, FMSI optionally constructs the kLCP array (a binary version of the truncated LCP array) 
to allow *O(1)* time for a *k*-mer in streaming queries.

The memory consumption of the index is split as follows:
- $2.5$ bits per superstring character to store the BWT ($2.5 - 3.5$ bits per *k*-mer for typical genomes and pan-genomes).
- Typically $0 - 0.8$ bits per *k*-mer to store the SA-transformed mask.
- Optionally $1$ bit per superstring character to store the kLCP array.
- $1$ byte per character of the largest entry in the queried fasta file (which is typically negligible).

The construction of the index requires the full suffix array to be computed and requires up to 17 bytes per superstring character.


## How to test

To run the associated tests, simply run `make test`.

## Issues

Please use [GitHub issues](https://github.com/OndrejSladky/fmsi/issues).

## Changelog

See [Releases](https://github.com/OndrejSladky/fmsi/releases).


## Licence

[MIT](https://github.com/OndrejSladky/fmsi/blob/master/LICENSE.txt)


## Contact

[Ondrej Sladky](https://iuuk.mff.cuni.cz/~sladky/) \<ondra.sladky@gmail.com\>\

