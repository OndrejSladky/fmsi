#! /usr/bin/env python3
"""Answers of FMSI 2 against an independent reference that enumerates canonical k-mers.

The reference knows nothing of the index: from a masked superstring (S, M) it lists the set
K = {canonical(S[i:i+k]) : M[i] = 1} and answers every query from this set. Small masked superstrings are
generated deterministically (fixed seeds), with repeated, masked-out and palindromic k-mers, in three kinds
of masks:

  min-one    exactly one set occurrence per k-mer of K (in either orientation)
  max-one    every occurrence of every k-mer of K, in both orientations
  other      a random non-empty subset of the occurrences of every k-mer of K

Every index (presets --extra-small, --small, --fast; without and with the kLCP array; the compact file;
an index made again from `fmsi export`) is queried in every way (isolated k-mers and reads; unstreamed,
streamed with and without the kLCP array; one query at a time and in batches; one thread and four; the
layouts in memory), and every answer is compared with the reference:

  membership  the string of 0/1 of `fmsi query`, window by window: k-mers in either orientation, windows
              with a character other than ACGT (absent), sequences shorter than k, lower-case queries
  lookup      -1 exactly for the absent windows. The three kinds of hashes are checked apart:
              * canonical minimal perfect hash (one strand, min-one mask): a bijection K -> [0, |K|), the
                same value for a k-mer and its reverse complement;
              * hashes over set occurrences (one strand, any other mask): values in [0, set bits), equal
                for the two orientations of a k-mer only if the search finds the same occurrence; not minimal;
              * orientation-specific hashes (--fast): values in [0, number_of_hashes), injective on
                oriented k-mers for a min-one mask, a k-mer and its reverse complement apart.
  access      access(lookup(q)) is q or its reverse complement (for --fast: q itself); over all hashes,
              access lists the k-mers of the set occurrences; a value that is not a hash gives `-`
  mask type   what `fmsi info` reports never contradicts the mask (min-one is never claimed for a mask
              that is not)

With --legacy, FMSI 0.4 (`fmsi legacy` of this binary, or a binary of FMSI 0.4 with --legacy-binary) is
checked on the same inputs: its scan (any mask), its first-bit test -O (max-one masks) and its streaming -S.

Usage: ./reference_test.py [--fmsi ../fmsi] [--quick] [--legacy] [--legacy-binary PATH] [--seed N] [-v]
Exit status 1 if any answer differs.
"""

import argparse
import collections
import itertools
import os
import random
import shutil
import subprocess
import sys
import tempfile

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMPLEMENT)[::-1]


def canonical(s):
    r = rc(s)
    return s if s <= r else r


# 1. The reference ------------------------------------------------------------------------------------------

class Reference:
    """A masked superstring and what it represents, by explicit enumeration."""

    def __init__(self, superstring, mask, k):
        self.s, self.m, self.k = superstring, mask, k
        self.set_positions = [i for i in range(len(superstring) - k + 1) if mask[i]]
        self.kmers = {canonical(superstring[i:i + k]) for i in self.set_positions}
        occurrences = collections.Counter(canonical(superstring[i:i + k]) for i in self.set_positions)
        all_occurrences = collections.Counter(canonical(superstring[i:i + k])
                                              for i in range(len(superstring) - k + 1))
        self.min_one = all(count == 1 for count in occurrences.values())
        self.max_one = all(occurrences[kmer] == all_occurrences[kmer] for kmer in self.kmers)
        # The same per strand, a k-mer and its reverse complement taken as different strings: what an index
        # of one strand can compute (`fmsi info`: mask_type_source computed-one-strand).
        oriented = collections.Counter(superstring[i:i + k] for i in self.set_positions)
        all_oriented = collections.Counter(superstring[i:i + k] for i in range(len(superstring) - k + 1))
        self.min_one_strand = all(count == 1 for count in oriented.values())
        self.max_one_strand = all(oriented[kmer] == all_oriented[kmer] for kmer in oriented)

    def fasta(self):
        return ">ms\n" + "".join(c if bit else c.lower() for c, bit in zip(self.s, self.m)) + "\n"

    def window_present(self, window):
        window = window.upper()
        return all(c in "ACGT" for c in window) and canonical(window) in self.kmers

    def membership(self, sequence):
        k = self.k
        return "".join("1" if self.window_present(sequence[i:i + k]) else "0"
                       for i in range(len(sequence) - k + 1))


def random_reference(rng, k, length, mask_kind, density):
    """A random superstring with repeats (small alphabets of k-mers collide; segments are repeated, also
    reverse-complemented) and a mask of the given kind for a random subset of its k-mers."""
    pieces = []
    total = 0
    while total < length:
        if pieces and rng.random() < 0.3:
            piece = rng.choice(pieces)
            piece = piece[rng.randrange(len(piece)):][:rng.randrange(1, 3 * k + 2)] or piece
            if rng.random() < 0.5:
                piece = rc(piece)
        else:
            piece = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1, 4 * k + 2)))
        if k % 2 == 0 and rng.random() < 0.3:                 # a palindromic k-mer
            half = "".join(rng.choice("ACGT") for _ in range(k // 2))
            piece += half + rc(half)
        pieces.append(piece)
        total += len(piece)
    s = "".join(pieces)
    if len(s) < k:
        s += "".join(rng.choice("ACGT") for _ in range(k - len(s)))
    positions = collections.defaultdict(list)
    for i in range(len(s) - k + 1):
        positions[canonical(s[i:i + k])].append(i)
    chosen = [kmer for kmer in sorted(positions) if rng.random() < density] or [sorted(positions)[0]]
    mask = [False] * len(s)
    for kmer in chosen:
        occurrences = positions[kmer]
        if mask_kind == "min-one":
            picked = [rng.choice(occurrences)]
        elif mask_kind == "max-one":
            picked = occurrences
        else:
            picked = [i for i in occurrences if rng.random() < 0.5] or [rng.choice(occurrences)]
        for i in picked:
            mask[i] = True
    return Reference(s, mask, k)


def queries_for(rng, ref, exhaustive_limit):
    """(name, sequence) pairs: isolated k-mers (all 4^k of them for small k; else k-mers of the superstring in
    both orientations, mutated ones and random ones) and reads (substrings of the superstring across masked
    stretches, in both orientations, with substitutions, with N and other characters, in lower case, shorter
    than k, empty, and longer random ones)."""
    k, s = ref.k, ref.s
    kmers = []
    if 4 ** k <= exhaustive_limit:
        kmers = ["".join(p) for p in itertools.product("ACGT", repeat=k)]
    else:
        for i in range(len(s) - k + 1):
            kmers += [s[i:i + k], rc(s[i:i + k])]
        for _ in range(200):
            i = rng.randrange(len(s) - k + 1)
            j = rng.randrange(k)
            kmer = s[i:i + k]
            kmers.append(kmer[:j] + rng.choice("ACGT") + kmer[j + 1:])
            kmers.append("".join(rng.choice("ACGT") for _ in range(k)))
    reads = [s, rc(s), "", "A", s[:k - 1], "N" * k, "N" * (2 * k + 1), s[:k] + "N", "N" + s[:k]]
    for _ in range(60):
        i = rng.randrange(len(s))
        read = s[i:i + rng.randrange(1, 6 * k + 3)]
        if rng.random() < 0.5:
            read = rc(read)
        kind = rng.randrange(6)
        if kind == 1 and read:                                    # substitutions
            read = list(read)
            for _ in range(rng.randrange(1, 4)):
                read[rng.randrange(len(read))] = rng.choice("ACGT")
            read = "".join(read)
        elif kind == 2 and read:                                  # invalid characters
            read = list(read)
            for _ in range(rng.randrange(1, 4)):
                read[rng.randrange(len(read))] = rng.choice("NNNRYKM-*")
            read = "".join(read)
        elif kind == 3:                                           # lower case
            read = read.lower()
        elif kind == 4:                                           # two pieces of the superstring joined
            j = rng.randrange(len(s))
            read += rc(s[j:j + 3 * k])
        elif kind == 5:                                           # an invalid character close to the end
            read = read + "N" + read[:rng.randrange(0, k)]
        reads.append(read)
    reads += ["".join(rng.choice("ACGT") for _ in range(8 * k)) for _ in range(5)]
    named = [("k%d" % i, kmer) for i, kmer in enumerate(kmers)]
    named += [("r%d" % i, read) for i, read in enumerate(reads)]
    return named


# 2. Running FMSI ---------------------------------------------------------------------------------------------

class Failures:
    def __init__(self, verbose):
        self.count = collections.Counter()
        self.checks = collections.Counter()
        self.examples = []
        self.verbose = verbose

    def check(self, category, condition, message):
        self.checks[category] += 1
        if not condition:
            self.count[category] += 1
            if self.count[category] <= 4:
                self.examples.append("%s: %s" % (category, message))
            if self.verbose:
                print("FAIL %s: %s" % (category, message), file=sys.stderr)
        return condition


def run(command, stdin=None):
    result = subprocess.run(command, input=stdin, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def parse_rows(stdout):
    """header -> answer of each output line (`fmsi query`, `fmsi lookup`)."""
    rows = []
    for line in stdout.split("\n")[:-1]:
        header, _, answer = line.partition("\t")
        rows.append((header, answer))
    return rows


def info_of(fmsi, index):
    code, out, _ = run([fmsi, "info", index])
    return dict(line.split("\t", 1) for line in out.splitlines() if "\t" in line and not line.startswith("section"))


def check_membership(fails, category, context, ref, queries, stdout, code):
    if not fails.check(category, code == 0, "%s: exit status %d" % (context, code)):
        return
    rows = parse_rows(stdout)
    if not fails.check(category, len(rows) == len(queries),
                       "%s: %d output lines for %d queries" % (context, len(rows), len(queries))):
        return
    for (name, sequence), (header, answer) in zip(queries, rows):
        expected = ref.membership(sequence)
        if not fails.check(category, header == name and answer == expected,
                           "%s: query %s %r: expected %s, got %s %s" % (context, name, sequence, expected,
                                                                         header, answer)):
            return


def parse_hashes(answer):
    return [int(x) for x in answer.split(",")] if answer else []


def check_lookup(fails, context, ref, queries, stdout, code, strands, n_hashes, claims_min_one):
    """Checks the hashes of one run; returns {oriented window: hash} for the present windows."""
    category = "lookup"
    hashes = {}
    if not fails.check(category, code == 0, "%s: exit status %d" % (context, code)):
        return hashes
    rows = parse_rows(stdout)
    if not fails.check(category, len(rows) == len(queries), "%s: %d output lines for %d queries"
                       % (context, len(rows), len(queries))):
        return hashes
    k = ref.k
    for (name, sequence), (header, answer) in zip(queries, rows):
        expected = ref.membership(sequence)
        values = parse_hashes(answer)
        got = "".join("0" if value < 0 else "1" for value in values)
        if not fails.check(category, header == name and got == expected,
                           "%s: query %s %r: membership by hashes %s, expected %s" % (context, name, sequence,
                                                                                    got, expected)):
            return hashes
        for i, value in enumerate(values):
            if value < 0:
                continue
            window = sequence[i:i + k].upper()
            fails.check(category, value < n_hashes, "%s: hash %d of %s outside [0, %d)" % (context, value, window, n_hashes))
            if window in hashes and ref.min_one:
                # A function of the k-mer; except a palindrome (even k) on both strands, which has two set
                # occurrences and thus two hashes (documented: Section S1.6).
                if strands == 2 and window == rc(window):
                    fails.check("palindrome-two-hashes", True, "")
                else:
                    fails.check("lookup-deterministic", hashes[window] == value,
                                "%s: %s has hashes %d and %d" % (context, window, hashes[window], value))
            hashes.setdefault(window, value)
    return hashes


def check_hash_semantics(fails, context, ref, hashes, strands, n_hashes):
    """The three kinds of hashes, each by its own definition."""
    by_canonical = collections.defaultdict(set)
    for window, value in hashes.items():
        by_canonical[canonical(window)].add(value)
    if strands == 1 and ref.min_one:
        # Canonical minimal perfect hash: one value per canonical k-mer, a bijection onto [0, |K|).
        fails.check("hash-canonical-mphf", n_hashes == len(ref.kmers),
                    "%s: %d hashes for %d k-mers" % (context, n_hashes, len(ref.kmers)))
        fails.check("hash-canonical-mphf", all(len(values) == 1 for values in by_canonical.values()),
                    "%s: a k-mer and its reverse complement have different hashes" % context)
        used = [next(iter(values)) for values in by_canonical.values()]
        fails.check("hash-canonical-mphf", len(set(used)) == len(used), "%s: two k-mers share a hash" % context)
        fails.check("hash-canonical-mphf", all(0 <= value < len(ref.kmers) for value in used),
                    "%s: a hash outside [0, |K|)" % context)
    elif strands == 1:
        # Hashes over set occurrences: within [0, set bits); distinct k-mers never share one.
        fails.check("hash-set-occurrences", n_hashes == len(ref.set_positions),
                    "%s: %d hashes for %d set bits" % (context, n_hashes, len(ref.set_positions)))
        owner = {}
        for kmer, values in by_canonical.items():
            for value in values:
                fails.check("hash-set-occurrences", owner.setdefault(value, kmer) == kmer,
                            "%s: hash %d shared by %s and %s" % (context, value, owner.get(value), kmer))
    else:
        # Orientation-specific hashes of the index of both strands.
        owner = {}
        for window, value in hashes.items():
            fails.check("hash-oriented", owner.setdefault(value, window) == window,
                        "%s: hash %d shared by %s and %s" % (context, value, owner.get(value), window))
        if ref.min_one:
            palindromes = sum(1 for kmer in ref.kmers if kmer == rc(kmer))
            fails.check("hash-oriented", n_hashes == 2 * len(ref.kmers),
                        "%s: %d hashes for %d k-mers (%d palindromes)" % (context, n_hashes, len(ref.kmers),
                                                                          palindromes))


def check_access(fails, context, fmsi, index, ref, hashes, strands, n_hashes, batch):
    category = "access"
    values = list(range(n_hashes)) + [n_hashes, n_hashes + 7]
    code, out, _ = run([fmsi, "access", "-b", str(batch), index], "".join("%d\n" % v for v in values))
    if not fails.check(category, code == 0, "%s: exit status %d" % (context, code)):
        return
    lines = out.split("\n")[:-1]
    if not fails.check(category, len(lines) == len(values), "%s: %d lines for %d hashes" % (context, len(lines), len(values))):
        return
    fails.check(category, lines[n_hashes:] == ["-", "-"], "%s: values beyond the hashes give %r" % (context, lines[n_hashes:]))
    kmers = lines[:n_hashes]
    k = ref.k
    set_occurrences = collections.Counter(ref.s[i:i + k] for i in ref.set_positions)
    if strands == 1:
        # Access lists the k-mers of the set occurrences, as they occur in the superstring.
        fails.check(category, collections.Counter(kmers) == set_occurrences,
                    "%s: access over all hashes does not list the set occurrences" % context)
    else:
        both = set_occurrences + collections.Counter(rc(kmer) for kmer in set_occurrences.elements())
        fails.check(category, collections.Counter(kmers) == both,
                    "%s: access over all hashes does not list the set occurrences of both strands" % context)
    # Round trip: access(lookup(q)) is q (both strands) or q in either orientation (one strand).
    for window, value in hashes.items():
        if value >= n_hashes:
            continue
        back = kmers[value]
        good = back == window if strands == 2 else canonical(back) == canonical(window)
        if not fails.check("round-trip", good, "%s: access(lookup(%s) = %d) = %s" % (context, window, value, back)):
            break
    # Round trip the other way for a min-one mask: lookup(access(h)) = h.
    if ref.min_one:
        for value, kmer in enumerate(kmers):
            if strands == 2 and kmer == rc(kmer):
                continue                      # a palindrome on both strands has two hashes (Section S1.6)
            if kmer in hashes and not fails.check("round-trip", hashes[kmer] == value,
                                                  "%s: lookup(access(%d) = %s) = %d" % (context, value, kmer, hashes[kmer])):
                break


QUERY_MODES = [("unstreamed", []), ("streamed", ["-S"]), ("streamed-noklcp", ["-S", "-L"])]


def test_reference(fails, fmsi, ref, queries, directory, tag, quick, legacy_binaries, rng):
    fasta = os.path.join(directory, tag + ".ms.fa")
    write(fasta, ref.fasta())
    query_file = os.path.join(directory, tag + ".q.fa")
    write(query_file, "".join(">%s\n%s\n" % (name, sequence) for name, sequence in queries))
    k = ref.k

    builds = [(preset, klcp, extra) for preset in ("extra-small", "small", "fast") for klcp in (False, True)
              for extra in ([],)]
    builds += [("small", False, ["-z"]), ("fast", True, ["-z"]), ("fast", True, ["--lean"]), ("small", True, ["-Q"]),
               ("small", True, ["-t", "4"]), ("fast", True, ["-t", "4"])]
    if quick:
        builds = [b for b in builds if not b[2] or b[2] == ["-z"]]
    reference_files = {}
    for number, (preset, klcp, extra) in enumerate(builds):
        name = "%s.%d.%s%s" % (tag, number, preset, "+klcp" if klcp else "")
        prefix = os.path.join(directory, name)
        index = prefix + ".fmsi"
        context = "%s k=%d %s%s %s" % (tag, k, preset, "+klcp" if klcp else "", " ".join(extra))
        command = [fmsi, "index", "--" + preset, "-k", str(k), "-o", index] + (["-l"] if klcp else []) + extra + [fasta]
        code, _, err = run(command)
        if k == 1 and klcp and preset != "fast":
            # No kLCP array for k = 1 on one strand: refused with an error message, not a crash.
            fails.check("build", code == 1 and "ERROR" in err, "%s: exit status %d" % (context, code))
            continue
        if not fails.check("build", code == 0 and os.path.exists(index), "%s: %s" % (context, err.strip()[-300:])):
            continue
        fails.check("build", os.path.exists(index + ".klcp") == klcp, "%s: kLCP file" % context)
        info = info_of(fmsi, index)
        strands = int(info.get("strands", 0))
        fails.check("info", strands == (2 if preset == "fast" else 1) and int(info.get("k", -1)) == k,
                    "%s: strands %s, k %s" % (context, info.get("strands"), info.get("k")))
        mask_type = info.get("mask_type", "")
        source = info.get("mask_type_source", "")
        # The type computed by an index of both strands is exact. An index of one strand computes it per
        # strand (source computed-one-strand): a mask that sets a k-mer once in each orientation is then
        # reported min-one although its hashes are not a canonical minimal perfect hash; this known limit is
        # counted (mask-type-per-strand-only), and the hashes are checked by what the mask really is.
        exact = strands == 2
        single = ref.min_one if exact else ref.min_one_strand
        uniform = ref.max_one if exact else ref.max_one_strand
        expected = "min-one,max-one" if single and uniform else "min-one" if single else "max-one" if uniform else "other"
        fails.check("mask-type", mask_type == expected and source == ("computed" if exact else "computed-one-strand"),
                    "%s: type %r (%s), expected %r" % (context, mask_type, source, expected))
        if not exact and (("min-one" in mask_type and not ref.min_one) or ("max-one" in mask_type and not ref.max_one)):
            fails.check("mask-type-per-strand-only", True, "")
        # Options that must not change the index (threads, the lean construction, QSufSort): equal files.
        key = (preset, klcp)
        if not extra:
            reference_files[key] = index
        elif extra != ["-z"] and key in reference_files:
            same = open(index, "rb").read() == open(reference_files[key], "rb").read()
            if klcp:
                same = same and open(index + ".klcp", "rb").read() == open(reference_files[key] + ".klcp", "rb").read()
            # The header holds no time stamp; the files of the same index are equal byte by byte.
            fails.check("same-file", same, "%s: file differs from the default construction" % context)

        # Membership: every way of querying.
        layouts = [[]] if quick else [[], ["-c", "64", "-m", "plain"], ["-c", "128", "-m", "packed"],
                                      ["-c", "256", "-m", "rrr"], ["-p", "0"]]
        for (mode, flags), batch, threads, layout in itertools.product(QUERY_MODES, (1, 32), (1, 4), layouts):
            if layout and (batch == 1 or threads == 4):
                continue
            code, out, _ = run([fmsi, "query"] + flags + layout + ["-b", str(batch), "-t", str(threads), "-q",
                                                                    query_file, index])
            category = "membership-" + mode
            check_membership(fails, category, "%s %s -b %d -t %d %s" % (context, mode, batch, threads,
                                                                        " ".join(layout)), ref, queries, out, code)
        # A small block of sequences per thread exercises the hand-over between threads.
        code, out, _ = run([fmsi, "query", "-S", "-t", "3", "-B", "64", "-q", query_file, index])
        check_membership(fails, "membership-threads", context + " -S -t 3 -B 64", ref, queries, out, code)

        # Dictionary: lookup in every way, then access.
        code, out, err = run([fmsi, "access", index], "0\n")
        n_hashes = None
        for cell in err.split():
            if cell.startswith("number_of_hashes="):
                n_hashes = int(cell.split("=")[1])
        if not fails.check("access", n_hashes is not None, "%s: no number_of_hashes" % context):
            continue
        first = None
        for (mode, flags), batch, threads in itertools.product(QUERY_MODES, (1, 32), (1, 4)):
            code, out, err = run([fmsi, "lookup"] + flags + ["-b", str(batch), "-t", str(threads), "-q", query_file, index])
            this = "%s lookup %s -b %d -t %d" % (context, mode, batch, threads)
            # Hashes that are not known to be minimal are announced: a mask that is not min-one, and a mask
            # found min-one on one strand only (it may set a k-mer once in each orientation).
            if first is None:
                announced = "Note: the mask" in err
                fails.check("lookup-note", announced == (mask_type not in ("min-one", "min-one,max-one")
                                                         or source == "computed-one-strand"),
                            "%s: note on the hashes %s" % (this, "printed" if announced else "missing"))
            hashes = check_lookup(fails, this, ref, queries, out, code, strands, n_hashes, "min-one" in mask_type)
            check_hash_semantics(fails, this, ref, hashes, strands, n_hashes)
            if first is None:
                first = hashes
            elif ref.min_one:
                # Hashes that are a function of the k-mer do not depend on how it is searched.
                fails.check("lookup-same-across-modes", hashes == first, "%s: hashes differ from the unstreamed ones" % this)
        for batch in (1, 32):
            check_access(fails, "%s access -b %d" % (context, batch), fmsi, index, ref, first or {}, strands, n_hashes, batch)

        # Saving and loading: the masked superstring comes back, and its index is the same file.
        if not extra:
            code, out, _ = run([fmsi, "export", index])
            fails.check("export", code == 0 and out.split("\n")[1:2] == [ref.fasta().split("\n")[1]],
                        "%s: export differs from the masked superstring" % context)

    # FMSI 0.4 on the same inputs: its scan (any mask), -O (max-one masks only), -S.
    for label, legacy in legacy_binaries:
        copy = os.path.join(directory, "%s.%s.fa" % (tag, label))
        shutil.copy(fasta, copy)
        code, _, err = run(legacy + ["index", "-k", str(k), copy])
        if not fails.check("legacy-build", code == 0, "%s %s: %s" % (tag, label, err.strip()[-200:])):
            continue
        plain = [(name, sequence) for name, sequence in queries
                 if sequence and all(c in "ACGTacgt" for c in sequence) and len(sequence) >= k]
        plain_file = os.path.join(directory, "%s.%s.q.fa" % (tag, label))
        write(plain_file, "".join(">%s\n%s\n" % (name, sequence) for name, sequence in plain))
        variants = [("legacy-scan", []), ("legacy-scan-streamed", ["-S"])]
        if ref.max_one:
            variants += [("legacy-first-bit", ["-O"]), ("legacy-first-bit-streamed", ["-O", "-S"])]
        for category, flags in variants:
            code, out, _ = run(legacy + ["query", "-k", str(k)] + flags + ["-q", plain_file, copy])
            check_membership(fails, category, "%s %s k=%d %s" % (tag, label, k, " ".join(flags)), ref, plain, out, code)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    parser.add_argument("--fmsi", default=os.path.join(here, "..", "fmsi"))
    parser.add_argument("--quick", action="store_true", help="fewer inputs and options (about a minute)")
    parser.add_argument("--legacy", action="store_true", help="also check `fmsi legacy` (the engine of FMSI 0.4)")
    parser.add_argument("--legacy-binary", help="also check a binary of FMSI 0.4")
    parser.add_argument("--seed", type=int, default=20261116)
    parser.add_argument("--keep", help="directory to keep the files in")
    parser.add_argument("--only", help="test only the inputs whose name contains this (e.g. r0.k31.other)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    fmsi = os.path.abspath(args.fmsi)
    legacy_binaries = []
    if args.legacy:
        legacy_binaries.append(("legacy", [fmsi, "legacy"]))
    if args.legacy_binary:
        legacy_binaries.append(("fmsi04", [os.path.abspath(args.legacy_binary)]))
    fails = Failures(args.verbose)
    rng = random.Random(args.seed)
    # (k, length of the superstring, share of the k-mers in the set); small k give exhaustive queries and
    # many repeats, even k palindromes, 31 and 32 the word boundary.
    shapes = [(1, 12, 0.7), (2, 40, 0.6), (3, 80, 0.5), (4, 200, 0.5), (5, 400, 0.6), (6, 600, 0.7),
              (8, 1500, 0.8), (11, 2500, 0.9), (13, 3000, 0.5), (16, 3000, 0.9), (21, 4000, 0.9),
              (31, 5000, 0.95), (32, 5000, 0.9), (31, 300, 0.3)]
    rounds = 1 if args.quick else 2
    if args.quick:
        shapes = [shape for shape in shapes if shape[0] in (2, 3, 4, 5, 8, 16, 31, 32)]
    directory = args.keep or tempfile.mkdtemp(prefix="fmsi-reference-test.")
    os.makedirs(directory, exist_ok=True)
    inputs = 0
    for round_ in range(rounds):
        for k, length, density in shapes:
            for mask_kind in ("min-one", "max-one", "other"):
                ref = random_reference(rng, k, length, mask_kind, density)
                queries = queries_for(rng, ref, 4 ** 6)
                tag = "r%d.k%d.%s" % (round_, k, mask_kind)
                if args.only and args.only not in tag:
                    continue                  # generated all the same, so that the other inputs do not change
                before = sum(fails.count.values())
                test_reference(fails, fmsi, ref, queries, directory, tag, args.quick, legacy_binaries, rng)
                inputs += 1
                status = "ok" if sum(fails.count.values()) == before else "FAILED"
                print("%-22s superstring %5d  k-mers %5d  set bits %5d  queries %5d  %s"
                      % (tag, len(ref.s), len(ref.kmers), len(ref.set_positions), len(queries), status), flush=True)
    if not args.keep:
        shutil.rmtree(directory)
    print("\n%-28s %9s %7s" % ("category", "checks", "failed"))
    for category in sorted(fails.checks):
        print("%-28s %9d %7d" % (category, fails.checks[category], fails.count[category]))
    failed = sum(fails.count.values())
    print("\n%d masked superstrings, %d checks, %d failed" % (inputs, sum(fails.checks.values()), failed))
    for example in fails.examples:
        print("  " + example)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
