#!/usr/bin/env python3

import argparse
from typing import List
import sys




def split_superstring(superstring: str, k: int) -> List[str]:
    """
    Split the given superstring into strings where at every position (except for the last k-1) begins a k-mer.
    """
    ret: List[str] = []
    next_start = 0
    was_present = False
    for i, c in enumerate(superstring):
        if c.islower():
            if was_present:
                ret.append(superstring[next_start: i + (k - 1)].upper())
            next_start = i + 1
            was_present = False
        else:
            was_present = True
    return ret


def get_k(superstring: str) -> int:
    """
    Determine the k given a masked superstring.
    """
    ret = 1
    while ret <= len(superstring):
        # If the k-th position from the end is the first with capital letter,
        # it means there starts the last k-mer of size k.
        if superstring[-ret].isupper():
            return ret
        ret += 1
    raise ValueError("Given superstring contains no valid k-mer.")


def main():
    """
    Read the superstring from standard input and print the rSPSS representation.
    """
    parser = argparse.ArgumentParser(description="Convert a masked superstring to rSPSS representation it induces.",
                                     usage="\nThe masked superstring is split at every position with 0 in the mask,\n"
                                     "which results in matchtigs (rSPSS) containing the same set of k-mers\n"
                                     "as was carried by the masked superstring (i.e. those k-mers masked with 1s)\n"
                                     "furthermore, if the mask minimized number of ones, the result are simplitigs (SPSS)."
                                     )
    parser.add_argument("infile", nargs="?", type=argparse.FileType("r"), default=sys.stdin)
    args = parser.parse_args()

    if not args.infile.readline().startswith(">"):
        print("ERROR: Fasta file should start with '>'.", file=sys.stderr)
        exit(1)

    superstring = args.infile.readline().strip()
    k = get_k(superstring)
    split = split_superstring(superstring, k)

    # Print the individual sequences.
    for i, s in enumerate(split):
        print(f">{i + 1}")
        print(s)


if __name__ == "__main__":
    main()
