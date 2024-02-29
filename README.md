# FMSI ($f$-Masked Superstring Index)

FMSI provides an experimental implementation of FMS-index - a BWT based indexing tool for
$f$-masked superstrings, which enables to ask membership queries and perform set operations on the
$k$-mer sets represented by the masked superstrings.

The implementation is based on the [SDSL library](https://github.com/simongog/sdsl-lite).


It is provided under the MIT license (see LICENSE file).

## How to cite

Coming soon!

## How to install

First clone the repo and its dependencies:

```
git clone --recursive git@github.com:OndrejSladky/fmsi.git
```

Compile the program by running `make`.

## How to run

To run the tool, run `./fmsi [command]`

The recognized commands are:

- `index` Creates a BWT based index of the given masked superstring.
- `query` Queries a $k$-mer against an index.
- `merge` Merge several indexes.
- `compact` Compacts an index.
- `export` Export the underlying masked superstring.
- `clean` Cleans the files stored for index.
- `-v`    Prints the version of the program.

### Index

Index (`./fmsi index`) recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta file with masked superstring to be indexed. This is a required argument.
- `-k value_of_k`    - The size of one k-mer. If not provided k is computed from the masked superstring under the assumption that the last run of zeros has length k-1.

For example: `./fmsi index -p spneumoniae.fa -k 13` 

### Query

Query (`./fmsi query`) returns whether the provided $k$-mer is in the masked superstring or not.
Note that `./fmsi index` must be run on the provided fasta file beforehand.

It recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta file from which the index was created. Required.
- `-q path_to_queries` - The path to the file with $k$-mer to be queried (each on a separate line). Required.
- `-k value_of_k`    - The size of one k-mer. Required.
- `-f function`      - A function to determine whether a $k$-mer
is represented based on the number of set and unset occurrences.
The recognized functions are following:
  - `or`  - Consider $k$-mer represented when any of its occurrences is set. This is the default function.
  - `all` - Assume that all occurrences are either set or unset and determine the presence by arbitrary occurrence.
  - `and` - Consider $k$-mer represented when all its occurrences are set.
  - `xor` - Consider $k$-mer represented when an odd number of occurrences is set.
  - `X-Y` (where X and Y can be any integers) - Consider $k$-mer represented when its number of set occurrences is between X and Y (inclusive).

For example: `./fmsi query -p spneumoniae.fa -f xor ACGT`

### Merge

Merge (`./fmsi merge`) creates a single index from several indexes representing the concatenation of the
corresponding masked superstrings.
Note that `./fmsi index` must be run on the provided fasta files beforehand.

It recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta file from which the index was created. Can be provided multiple times. Required at least twice.
- `-r path_to_result` - The path to the file where the result will be stored. Required.
- `-k value_of_k`    - The size of one k-mer.

For example: `./fmsi merge -p spneumoniae1.fa -p spneumoniae2.fa -r spneumoniae.fa -k 13`

### Compact

Compact (`./fmsi compacts`) compacts the index by removing the redundant information by computing a new or-masked superstring.
Especially note that the resulting superstring is or-masked superstring and not f-masked superstring.

It recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta file from which the index was created. This is a required argument.
- `-k value_of_k`    - The size of one k-mer. Required.
- `-f function`      - A function to determine whether a $k$-mer is represented based on the number of set and unset occurrences.
  The recognized functions are following:
  - `or`  - Consider $k$-mer represented when any of its occurrences is set. This is the default function.
  - `all` - Assume that all occurrences are either set or unset and determine the presence by arbitrary occurrence.
  - `and` - Consider $k$-mer represented when all its occurrences are set.
  - `xor` - Consider $k$-mer represented when an odd number of occurrences is set.
  - `X-Y` (where X and Y can be any integers) - Consider $k$-mer represented when its number of set occurrences is between X and Y (inclusive).
- `-s`                - Only print the masked superstring and do not modify the index.
- `-l`                - Use the local algorithm acting directly on the index.
- `-d value_of_d_max` - The maximum extension length. Relevant only for the local algorithm. Default is 5.

For example: `./fmsi compact -p spneumoniae.fa -k 13 -f xor -s`

### Export

Export (`./fmsi export`) recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta from which the index was created. This is a required argument.


### Clean

Clean (`./fmsi clean`) recognizes the following arguments:

- `-p path_to_fasta` - The path to the fasta from which the index was created. This is a required argument.


## How to run unittests

To run the associated tests, simply run `make test`.

