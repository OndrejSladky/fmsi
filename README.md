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
- `clean` Cleans the files stored for index.
- `-v`    Prints the version of the program.

### Index

Index (`./ms-index index`) recognizes the following arguments:

- `-k value_of_k` - the size of one k-mer. If not provided k is computed from the masked superstring under the assumption that the last run of zeros has length k-1.
- `-l value_of_l` - the size of l for which a mask is computed from the mask for k. l should not be greater than k. This argument can be provided multiple times.

Additionally, the last argument is the path to the fasta file containing the masked superstring
to be indexed.

For example: `./ms-index index -k 13 spneumoniae.fa` 

### Query

Index (`./ms-index index`) return whether the provided $k$-mer is in the masked superstring or not.
Note that `./ms-index index` must be run on the provided fasta file beforehand.

It recognizes the optional parameter `f` which is the function to determine whether a $k$-mer
is represented based on the number of set and unset occurrences.
The recognized functions are following:

- `default` -- Assume that all occurrence are either set or unset and determined the presence by arbitrary occurrence.
- `or` -- Consider $k$-mer represented when any of its occurence is set.
- `and` -- Consider $k$-mer represented when all its occurrences are set.
- `or` -- Consider $k$-mer represented an odd number of occurences is set.
- `X-Y` (where X and Y can be any integers)-- Consider $k$-mer represented when its number of set occurrences is between X and Y (inclusive).

It then takes two positional arguments, the fasta file and the queried $k$-mer.

For example: `./ms-index query spneumoniae.fa ACGT`

### Clean

Clean (`./ms-index clean`) recognizes one positional argument - the path to the fasta file.


## How to test

To run the associated tests, simply run `make test`.

