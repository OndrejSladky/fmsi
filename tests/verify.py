"""Verification script to test that FMSI returns the correct results for a possibly specified query file."""
import random
import subprocess
import argparse

fmsi_path = "../fmsi"

def reverse_complement(s):
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[c] for c in reversed(s))

def lmbda_uni(superstring, mask, kmer):
    ret = []
    k = len(kmer)
    for i in range(len(superstring)):
        if superstring[i:i+k] == kmer:
            ret.append(mask[i])
    return ret

def lmbda(superstring, mask, kmer):
    return lmbda_uni(superstring, mask, kmer) + lmbda_uni(superstring, mask, reverse_complement(kmer))

def f_or(occurrences):
    return any(occurrences)

def separate_mask_and_superstring(masked_superstring):
    mask = []
    superstring = []
    for c in masked_superstring:
        mask.append(c.isupper())
        superstring.append(c.upper())
    return "".join(superstring), mask

def read_masked_superstring(file):
    with open(file, 'r') as f:
        _ = f.readline().strip()
        masked_superstring = f.readline().strip()
    return masked_superstring

def run_fmsi_index(file):
    subprocess.run([fmsi_path, 'index', '-p', file])

def create_fmsi_process(file, k):
    return subprocess.Popen([fmsi_path, 'query', '-p', file, '-q', '-', '-k', str(k), '-F'], stdout=subprocess.PIPE, stdin=subprocess.PIPE)

def assert_correct_results(process, superstring, mask, kmers, k: int):
    for i, kmer in enumerate(kmers):
        if len(kmer) != k:
            print(f"Skipping kmer {kmer} of length {len(kmer)}")
            continue
        process.stdin.write(f"{kmer}\n".encode())
        process.stdin.flush()
        output = process.stdout.readline().decode().strip()
        assert output in ["FOUND", "NOT FOUND"], f"Unexpected output: {output} for kmer {kmer}"
        present_fmsi = output == "FOUND"
        present = f_or(lmbda(superstring, mask, kmer))
        assert present_fmsi == present, f"Expected {present}, got {present_fmsi} for kmer {kmer}"
        print(".", end="")
        if i % 100 == 99:
            print(f" {i+1}/{len(kmers)}")

def generate_random_kmers(k, superstring, num_queries):
    threshold = 0.5 * (len(superstring) > k)
    kmers = []
    for _ in range(num_queries):
        if random.random() < threshold:
            start = random.randint(0, len(superstring) - k)
            kmers.append(superstring[start:start+k])
        else:
            kmers.append("".join(random.choice("ACGT") for _ in range(k)))
    return kmers

def main():
    print("Started testing")
    parser = argparse.ArgumentParser("check whether FMSI gives correct answers for a given dataset and possible a query file (otherwise positive and negative queries are generated randomly)")
    parser.add_argument("path", help="path to the fasta file on which ./kmers is verified")
    parser.add_argument("--k", type=int, help="the value of k; required")
    parser.add_argument("--query_path", help="the path to the query file; if not specified, random queries are generated")
    parser.add_argument("--num_queries", type=int, help="the number of queries to generate if no query file is specified", default=1000)
    args = parser.parse_args()

    run_fmsi_index(args.path)
    superstring, mask = separate_mask_and_superstring(read_masked_superstring(args.path))
    if args.query_path:
        with open(args.query_path, 'r') as f:
            kmers = f.readlines()
    else:
        kmers = generate_random_kmers(args.k, superstring, args.num_queries)
    process = create_fmsi_process(args.path, args.k)
    assert_correct_results(process, superstring, mask, kmers, args.k)
    print("")

random.seed(42)
main()
