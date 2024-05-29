# FMSI ($f$-Masked Superstring Index)

FMSI is a memory-efficient tool for querying $k$-mer sets and performing set operations on them.
Internally it implements the FMS-index - a BWT based index for $f$-masked superstrings.

Note that currently the tool is still in development and especially the set operations are quite time-consuming.

It is provided under the MIT license (see LICENSE file).

## How to cite

If you use FMSI in your research, please cite the following.

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Function-Assigned Masked Superstrings as a Versatile and Compact Data Type for *k*-Mer Sets.
> *bioRxiv* 2024.03.06.583483, 2024. [https://doi.org/10.1101/2024.03.06.583483](https://doi.org/10.1101/2024.03.06.583483)

```
@article {sladky2024-f-masked-superstrings,
	author = {Ond{\v r}ej Sladk{\'y} and Pavel Vesel{\'y} and Karel B{\v r}inda},
	title = {Function-Assigned Masked Superstrings as a Versatile and Compact Data Type for 𝑘-Mer Sets},
	elocation-id = {2024.03.06.583483},
	year = {2024},
	doi = {10.1101/2024.03.06.583483},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/03/11/2024.03.06.583483},
	eprint = {https://www.biorxiv.org/content/early/2024/03/11/2024.03.06.583483.full.pdf},
	journal = {bioRxiv}
}

```

> Ondřej Sladký, Pavel Veselý, and Karel Břinda: Masked superstrings as a unified framework for textual *k*-mer set representations. *bioRxiv* 2023.02.01.526717, 2023.
[https://doi.org/10.1101/2023.02.01.526717](https://doi.org/10.1101/2023.02.01.526717)

```
@article{sladky2023-masked-superstrings,
  title   = { Masked superstrings as a unified framework for textual $k$-mer set representations },
  author  = { Sladk{\'y}, Ond{\v r}ej and Vesel{\'y}, Pavel and B{\v r}inda, Karel },
  journal = { bioRxiv },
  volume  = { 2023.02.01.526717 },
  year    = { 2023 },
  doi     = { 10.1101/2023.02.01.526717 }
}
```

## How to install

First clone the repo and its dependencies:

```
git clone --recursive git@github.com:OndrejSladky/fmsi.git
```

Compile the program by running `make`.

## How to use

### Basic usage

If you wish to perform set operations or answer membership queries without the need to understand the
details of the $f$-masked superstring framework, FMSI can manage the details for you.
This can be done in a few simple steps:
1. **Compute a masked superstring.**
   - This can be done by [KmerCamel🐫](https://github.com/OndrejSladky/kmercamel); simply run `kmercamel -c -k 31 -p fasta.fa -o ms.fa` (with appropriate values for `-k` and `-p`).
   - If you obtained the masked superstring in a different way, make sure it minimizes the number of ones; if you're unsure, you can use `kmercamel optimize -c -a zeros -k 31 -p ms_more_ones.fa -o ms.fa`. No need to optimize superstrings directly computed by KmerCamel🐫.
2. **Index the masked superstring** with `./fmsi index -p ms.fa`.
3. **Perform the set operations** you wish with `./fmsi name_of_the_operation -k 31 -p ms1.fa -p ms2.fa -r ms.fa`. Possible operations (use the names instead of `name_of_the_operation`) are:
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

### Advanced usage

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

## How to run unittests

To run the associated tests, simply run `make test`.

## Contact

Ondřej Sladký - `sladky@iuuk.mff.cuni.cz`
