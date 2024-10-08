#!/usr/bin/env python3
import subprocess
import sys
import os
import argparse


def verify_instance(fasta_path: str, k: int, algorithm: str, complements: bool, large: bool, masked_superstring: str) -> bool:
    """
    Check if running superstring algorithm on given fasta file produces the same set of k-mers as the original one.
    """
    with open("./bin/kmercamel.txt", "w") as k_mers:
        args = ["./kmercamel-large" if large else "./kmercamel"] + (["optimize"] if masked_superstring else []) + ["-p", (masked_superstring if masked_superstring != "" else fasta_path), "-k", f"{k}", "-a", algorithm]
        if complements:
            args.append("-c")
        subprocess.run(args, stdout=k_mers)
    with open("./bin/kmercamel.txt", "r") as k_mers:
        with open("./bin/converted.fa", "w") as converted:
            subprocess.run(["./convert_superstring.py"], stdin=k_mers, stdout=converted)
    # in result; in original sequence; in result without complements; in original without complements; in merged file
    stats = [{}, {}, {}]
    runs = [
        (0, "./bin/converted.fa", "converted", complements),
        (1, fasta_path, "original", complements),
    ]
    for i, path, result, pass_complements in runs:
        args = ["jellyfish", "count", "-m", f"{k}", "-s", "100M", "-o", f"./bin/{result}.jf", path]
        if pass_complements:
            args.insert(2, "-C")
        subprocess.run(args)
        with open(f"./bin/{result}_stats.txt", "w") as f:
            subprocess.run(["jellyfish", "stats", f"./bin/{result}.jf"], stdout=f)
        with open(f"./bin/{result}_stats.txt", "r") as f:
            for _ in range(4):
                key, value = f.readline().split()
                stats[i][key] = value
    # Count k-mers on merged file.
    subprocess.run(["jellyfish", "merge", "-o", f"./bin/merged.jf", "./bin/converted.jf", "./bin/original.jf"])
    with open(f"./bin/merged_stats.txt", "w") as f:
        subprocess.run(["jellyfish", "stats", f"./bin/merged.jf"], stdout=f)
    with open(f"./bin/merged_stats.txt", "r") as f:
        for _ in range(4):
            key, value = f.readline().split()
            stats[2][key] = value
    distinct_key = "Distinct:"
    total_key = "Total:"
    if stats[0][distinct_key] != stats[1][distinct_key] or stats[0][distinct_key] != stats[2][distinct_key]:
        print("F")
        print(f"Failed: k={k}: expected orginal_distinct_count={stats[1][distinct_key]}, result_distinct_count={stats[0][distinct_key]} and merged_distinct_count={stats[2][distinct_key]} to be equal.")
        return False
    elif complements and stats[0][distinct_key] != stats[0][total_key] and masked_superstring == "":
        print("W")
        print(f"Warning: k={k}: number of masked k-mers={stats[0][total_key]} is not minimal possible (minimum is {stats[0][distinct_key]}).")
    else:
        print(".", end="")
        sys.stdout.flush()
    return True


def main():
    # Initialize.
    if not os.path.exists("bin"):
        os.makedirs("bin")

    parser = argparse.ArgumentParser("check whether ./kmercamel(-large) (optimize) outputs a superstring which contains the same set of k-mers"
                                     "as the original sequence")
    parser.add_argument("path", help="path to the fasta file on which ./kmers is verified")
    parser.add_argument("--quick", help="if set do not check for full range of k", action="store_true")
    parser.add_argument("--superstring_path", help="the path to the masked superstring to test masks")
    parser.add_argument("--k", help="the value of k for mask verification")
    args = parser.parse_args()

    success = True
    if args.superstring_path is None:
        # Do the tests on superstring algoritms.
        for a in ["global", "local", "globalAC", "localAC", "streaming"]:
            print(f"Testing {a}:")
            for complements in [True, False]:
                for large in [True, False]:
                    for k in ( ([5, 8, 12] if a not in ["local", "global"] else [5, 8, 12, 17, 31, 43, 51, 63]) if args.quick else range(2, 64)):
                        if not large and k >= 32: continue
                        success &= verify_instance(args.path, k, a, complements, large, "")
                    print("")
    else:
        k = int(args.k)
        # Do the tests on mask optimization.
        for a in ["runs", "ones", "zeros", "runsapprox"]:
            print(f"Testing {a}:")
            for large in [True, False]:
                if not large and k >= 32: continue
                success &= verify_instance(args.path, k, a, True, large, args.superstring_path)
            print("")

    # Print status.
    if not success:
        print("Tests failed")
        exit(1)
    print("OK")

main()
