# MS Index

MS Index provides an efficient implementation of a BWT based indexing tool for
[masked superstrings](https://doi.org/10.1101/2023.02.01.526717).

## How to install

First clone the repo and its dependencies:

```
git clone --recursive https://github.com/OndrejSladky/kmercamel
```

Compile the program by running `make`.

## How to run

To run the tool, run `./ms-index [command]`

The recognized commands are:

- `index` Creates a BWT based index of the given masked superstring.
- `query` Queries a $k$-mer against an index.
- `-v`    Prints the version of the program.

### Index

Index (`./ms-index index`) recognizes the following arguments:

- `-k value_of_k` - the size of one k-mer.

Additionally, the last argument is the path to the fasta file containing the masked superstring
to be indexed.

For example: `./ms-index index -k 13 spneumoniae.fa` 

### Query

Index (`./ms-index index`) return whether the provided $k$-mer is in the masked superstring or not.
Note that `./ms-index index` must be run on the provided fasta file beforehand.

Index takes two positional arguments, the fasta file and the queried $k$-mer.

For example: `./ms-index query spneumoniae.fa ACGT`

## How to test

To run the associated tests, simply run `make test`.

