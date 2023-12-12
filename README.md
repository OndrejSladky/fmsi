# MS Index

MS Index provides an efficient implementation of a BWT based indexing tool for
[masked superstrings](https://doi.org/10.1101/2023.02.01.526717).

## How to install

First clone the repo and its dependencies:

```
git clone --recursive git@github.com:OndrejSladky/ms-index.git
```

Compile the program by running `make`.

## How to run

To run the tool, run `./ms-index [command]`

The recognized commands are:

- `index` Creates a BWT based index of the given masked superstring.
- `query` Queries a $k$-mer against an index.
- `clean` Cleans the files stored for index.
- `-v`    Prints the version of the program.

### Index

Index (`./ms-index index`) recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta file with masked superstring to be indexed. This is a required argument.
- `-k value_of_k`    - The size of one k-mer. If not provided k is computed from the masked superstring under the assumption that the last run of zeros has length k-1.
- `-l value_of_l`    - The size of l for which a mask is computed from the mask for k. l should not be greater than k. This argument can be provided multiple times.

For example: `./ms-index index -p spneumoniae.fa -k 13 -l 12 -l 11` 

### Query

Query (`./ms-index query`) returns whether the provided $k$-mer is in the masked superstring or not.
Note that `./ms-index index` must be run on the provided fasta file beforehand.

It recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta file from which the index was created.
- `-f function`      - A function to determine whether a $k$-mer
is represented based on the number of set and unset occurrences.
The recognized functions are following:
  - `or`  - Consider $k$-mer represented when any of its occurrences is set. This is the default function.
  - `all` - Assume that all occurrences are either set or unset and determine the presence by arbitrary occurrence.
  - `and` - Consider $k$-mer represented when all its occurrences are set.
  - `xor` - Consider $k$-mer represented when an odd number of occurrences is set.
  - `X-Y` (where X and Y can be any integers) - Consider $k$-mer represented when its number of set occurrences is between X and Y (inclusive).
It then takes one positional argument - the queried $k$-mer.

For example: `./ms-index query -p spneumoniae.fa -f xor ACGT`

### Clean

Clean (`./ms-index clean`) recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta from which the index was created. This is a required argument.


## How to test

To run the associated tests, simply run `make test`.

