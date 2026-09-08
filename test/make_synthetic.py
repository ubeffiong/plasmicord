#!/usr/bin/env python3
"""
make_synthetic.py -- generate a tiny synthetic outbreak set so the WHOLE
framework runs offline (no mob_recon, no Mash).

It writes:
  <outdir>/plasmids/<isolate>__<name>.fasta   one record per plasmid
  <outdir>/plasmid_index.tsv                   plasmid_id, isolate_id, length
  <outdir>/metadata.tsv                        isolate_id, chromosomal_cluster, location

Design (so results are predictable):
  Two "backbone" plasmids A and B (random DNA).
  iso1 (CC1): A
  iso2 (CC1): A, B
  iso3 (CC2): A'   (A with ~1% mutations -> same plasmid unit as A)
  iso4 (CC2): B'   (B with ~1% mutations -> same plasmid unit as B)
  iso5 (CC3): none (plasmid-free)
Expected: PU for A shared across CC1 & CC2 (cross-link), PU for B likewise.
"""
import argparse
import os
import random

BASES = "ACGT"


def rand_seq(n, rng):
    return "".join(rng.choice(BASES) for _ in range(n))


def mutate(seq, rate, rng):
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < rate:
            s[i] = rng.choice([b for b in BASES if b != s[i]])
    return "".join(s)


def write_fasta(path, header, seq, width=70):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    pdir = os.path.join(args.outdir, "plasmids")
    os.makedirs(pdir, exist_ok=True)

    A = rand_seq(1200, rng)
    B = rand_seq(900, rng)
    A3 = mutate(A, 0.01, rng)
    B4 = mutate(B, 0.01, rng)

    # (isolate, plasmid_name, sequence)
    plasmids = [
        ("iso1", "pA", A),
        ("iso2", "pA", A),
        ("iso2", "pB", B),
        ("iso3", "pA", A3),
        ("iso4", "pB", B4),
        # iso5 has no plasmids
    ]

    index_path = os.path.join(args.outdir, "plasmid_index.tsv")
    with open(index_path, "w") as idx:
        idx.write("plasmid_id\tisolate_id\tlength\n")
        for iso, name, seq in plasmids:
            pid = f"{iso}__{name}"
            write_fasta(os.path.join(pdir, pid + ".fasta"), pid, seq)
            idx.write(f"{pid}\t{iso}\t{len(seq)}\n")

    meta_path = os.path.join(args.outdir, "metadata.tsv")
    with open(meta_path, "w") as md:
        md.write("isolate_id\tchromosomal_cluster\tlocation\n")
        for iso, cc in [("iso1", "CC1"), ("iso2", "CC1"),
                        ("iso3", "CC2"), ("iso4", "CC2"), ("iso5", "CC3")]:
            md.write(f"{iso}\t{cc}\tsite_{cc[-1]}\n")

    print(f"[make_synthetic] wrote plasmids to {pdir}")
    print(f"[make_synthetic] wrote {index_path} and {meta_path}")


if __name__ == "__main__":
    main()
