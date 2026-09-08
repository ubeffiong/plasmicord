#!/usr/bin/env python3
"""
kmer_distance.py -- dependency-free all-vs-all distance matrix over plasmid
FASTA files, for the offline demo and as a no-Mash fallback on SMALL inputs.

For each FASTA it builds the set of canonical k-mers (min of a k-mer and its
reverse complement), computes the Jaccard index J between two sets, and reports a
Mash-style distance:
        D = -1/k * ln( 2J / (1 + J) )      (0 when J = 1)
capped at 1.0 when J = 0.

This is EXACT Jaccard on full k-mer sets -- fine for a handful of small plasmids
(demo / teaching / tiny outbreak sets). For real datasets use Mash (sketching),
which is far faster and what the main pipeline calls. The output matrix format is
identical to mash_to_matrix.py, so downstream steps don't care which produced it.

Usage:
  kmer_distance.py --fastas f1.fasta f2.fasta ... --k 21 --out matrix.tsv
  kmer_distance.py --dir results/plasmids --k 21 --out matrix.tsv

Standard library only.
"""

import argparse
import glob
import math
import os
import sys

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta_concat(path):
    """Concatenate all sequence (a plasmid file may have >1 contig)."""
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                seq.append("N")  # never create artificial k-mers across contig boundaries
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def canonical_kmers(seq, k):
    kmers = set()
    if len(seq) < k:
        return kmers
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if set(km) - set("ACGT"):
            continue
        rc = revcomp(km)
        kmers.add(km if km <= rc else rc)
    return kmers


def mash_distance(j, k):
    if j <= 0:
        return 1.0
    if j >= 1:
        return 0.0
    return -1.0 / k * math.log(2 * j / (1 + j))


def item_id(path):
    b = os.path.basename(path)
    for e in (".fasta", ".fa", ".fna"):
        if b.endswith(e):
            return b[: -len(e)]
    return b


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--fastas", nargs="+")
    g.add_argument("--dir")
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.k < 1:
        ap.error("k must be positive")

    if args.dir:
        files = sorted(glob.glob(os.path.join(args.dir, "*.fasta")) +
                       glob.glob(os.path.join(args.dir, "*.fa")) +
                       glob.glob(os.path.join(args.dir, "*.fna")))
    else:
        files = args.fastas
    if not files:
        sys.exit("ERROR: no FASTA files found.")

    ids = [item_id(f) for f in files]
    if len(ids) != len(set(ids)):
        sys.exit("ERROR: duplicate plasmid filename identifiers.")
    ksets = [canonical_kmers(read_fasta_concat(f), args.k) for f in files]
    if any(not s for s in ksets):
        sys.exit("ERROR: each plasmid must contain at least one unambiguous k-mer.")
    n = len(files)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = ksets[i], ksets[j]
            union = len(a | b)
            inter = len(a & b)
            J = (inter / union) if union else 0.0
            d = mash_distance(J, args.k)
            mat[i][j] = mat[j][i] = d

    with open(args.out, "w") as out:
        out.write("item\t" + "\t".join(ids) + "\n")
        for i in range(n):
            out.write(ids[i] + "\t" + "\t".join(f"{mat[i][j]:.6f}" for j in range(n)) + "\n")

    sys.stderr.write(f"[kmer_distance] wrote {n}x{n} matrix (k={args.k}) -> {args.out}\n")


if __name__ == "__main__":
    main()
