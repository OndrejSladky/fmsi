# KmerCamel🐫 internals 

## Parsing FASTA files and storing k-mers

To parse FASTA files, we use the `kseq.h` library. To support large sequences, we use the version from [seqtk](https://github.com/lh3/seqtk/blob/master/kseq.h).
We represent *k*-mers as integers, where the size of the integer is selected depending on *k* without any need to recompile.
We use 64bit integers, 128bit integers from GCC and 256bit integers implemented in `uint256_t` folder.
To achieve this while keeping high performance, we use C++ templates and, where needed, C macros.
Efficient operations on *k*-mers are implemented in the `kmer.h` file.
*k*-mers are stored in a `khash.h` hash table. We modified the original version to internally use 64bit integers to support very large *k*-mer sets and also use Wang hash instead of the default one.
We implement wrapper operations over `khash.h` in `khash_utils.h` and the *k*-mer parser in `parser.h`.

## Global greedy

The global greedy algorithm in its core works in the unidirectional model as follows:
```python
def global_greedy(K):
    while len(K) > 1:
        a, b = two_most_overlapping(K)
        K = K + {merge(a, b)} - {a, b}
    return K.first()
```

In the bidirectional model, we additionally always merge the reverse complements of `a` and `b`, while never merging two reverse complementary strings.

To achieve high performance, the problem is split into two parts. In the first, we for each *k*-mer compute its successor
without actually merging them, and we obtain the masked superstring in the second part.
To quickly find the two most overlapping *k*-mers, we, starting from the largest overlap length,
create a map of prefixes to *k*-mers and then iterate over the suffixes.
Since this map is quite memory demanding, we store at each time only a part of the *k*-mers and repeat the process that many times (which can be turned off).
In order for this not to be as time-consuming, we first partially sort the *k*-mers using bucket sort.

The global greedy is implemented in the `global.h` file.

## Local greedy

The local greedy in its core works as follows:
```python
def local_greedy(K):
    result = ""
    while len(K) > 0:
        current = K.first()
        K = K - {current}
        while largest_overlap(K, current) >= k - d_max:
            a = most_overlapping(K, current)
            current = merge(a, current)
            K = K - {a}
        result = result + current
    return result
```

The *k*-mers with largest overlap are found simply by iterating over all possible extensions of length up to `d_max`.

The local greedy is implemented in the `local.h` file.

## Mask optimization

We implement three different mask optimization algorithms:
- Minimizing the number of 1s in the mask.
- Maximizing the number of 1s in the mask.
- Minimizing the number of runs of 1s in the mask.

Minimization/maximization of 1s is implemented by a simple two pass algorithm, where in the first pass *k*-mers
are loaded and in the second, the mask is masked at all positions or at the first position respectively.

Minimization of runs of 1s is done in three steps. First, the intervals of consecutive *k*-mers which can be masked to 1 are obtained.
Second, we set intervals with *k*-mers not appearing elsewhere to 1 and after all intervals with resolved *k*-mers to 0.
This step is important as it significantly reduces the problem size.
Finally, e pass the rest to a GLPK ILP solver to minimize the number of 1s in the mask.

The mask optimization is implemented in the `masks.h` file.

## Aho-Corasick-based versions and other experimental algorithms

All the experimental algorithms are in the folder `ac`.