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

def create_fmsi_process(file, k, optimize):
    return subprocess.Popen([fmsi_path, 'query', '-p', file, '-q', '-', '-k', str(k), '-F'] + (['-O'] if optimize else []), stdout=subprocess.PIPE, stdin=subprocess.PIPE)

def get_total_kmers(sequence, k):
    return len(sequence) - k + 1

def get_valid_kmers(sequence, k):
    return sum(all(c in "ACGT" for c in sequence[i:i+k]) for i in range(len(sequence) - k + 1))

def get_positive_kmers(sequence, superstring, mask, k):
    return sum(f_or(lmbda(superstring, mask, sequence[i:i+k])) for i in range(len(sequence) - k + 1))


def assert_correct_results(process, superstring, mask, kmers, k: int):
    results = []
    for i, kmer in enumerate(kmers):
        process.stdin.write(f">\n{kmer}\n".encode())
        results.append([get_total_kmers(kmer, k), get_valid_kmers(kmer, k), get_positive_kmers(kmer, superstring, mask, k)])
    process.stdin.close()
    index = 0
    for line in process.stdout:
        total_count, valid_count, positive_count = results[index]
        total_kmers, valid_kmers, found_kmers = map(int, line.decode().strip().split(","))
        assert total_count == total_kmers, f"Expected {total_kmers} kmers, got {total_count}"
        assert valid_count == valid_kmers, f"Expected {valid_kmers} valid kmers, got {valid_count}"
        assert positive_count == found_kmers, f"Expected {found_kmers} positive kmers, got {positive_count}"
        print(".", end="")
        if (index + 1) % 50 == 0:
            print()
        index += 1

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
    parser.add_argument("--streaming_length", type=int, help="the number of kmers in a single line (default 1)", default=1)
    parser.add_argument("--query_path", help="the path to the query file; if not specified, random queries are generated")
    parser.add_argument("--num_queries", type=int, help="the number of queries to generate if no query file is specified", default=1000)
    parser.add_argument("--optimize", type=bool, help="if queries should optimize for maximizing the number of ones in the mask", default=False)
    args = parser.parse_args()

    run_fmsi_index(args.path)
    superstring, mask = separate_mask_and_superstring(read_masked_superstring(args.path))
    if args.query_path:
        with open(args.query_path, 'r') as f:
            kmers = f.readlines()
    else:
        kmers = generate_random_kmers(args.k + args.streaming_length - 1, superstring, args.num_queries)
    process = create_fmsi_process(args.path, args.k, args.optimize)
    assert_correct_results(process, superstring, mask, kmers, args.k)
    print()
    print("OK")

random.seed(42)
main()
