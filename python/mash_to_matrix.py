#!/usr/bin/env python3
"""
mash_to_matrix.py -- convert `mash dist` output into the symmetric TSV matrix
that cluster_plasmids.py expects.

`mash dist ref.msh query.msh` prints one line per pair:
    <ref_path>  <query_path>  <distance>  <p_value>  <shared_hashes>

We map each path to an item id = basename with the FASTA extension stripped, and
build a full symmetric matrix (diagonal 0).

Usage:
  mash_to_matrix.py --dist mash_dist.tsv --out plasmid_matrix.tsv
Optional:
  --strip-ext ".fasta,.fa,.fna"   extensions to strip from filenames (default set)

Standard library only.
"""

import argparse
import math
import os
import sys


def item_id(path, exts):
    b = os.path.basename(path)
    for e in exts:
        if b.endswith(e):
            return b[: -len(e)]
    return b


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dist", required=True, help="mash dist output TSV")
    ap.add_argument("--out", required=True)
    ap.add_argument("--strip-ext", default=".fasta,.fa,.fna,.fasta.gz,.fna.gz")
    args = ap.parse_args()

    exts = [e for e in args.strip_ext.split(",") if e]

    pairs = {}
    items = []
    seen = set()
    with open(args.dist) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                sys.exit("ERROR: malformed Mash distance row.")
            a = item_id(f[0], exts)
            b = item_id(f[1], exts)
            d = float(f[2])
            if not math.isfinite(d) or not 0 <= d <= 1:
                sys.exit("ERROR: invalid Mash distance.")
            if (b, a) in pairs and abs(pairs[(b, a)] - d) > 1e-9:
                sys.exit("ERROR: asymmetric Mash distances.")
            pairs[(a, b)] = d
            for x in (a, b):
                if x not in seen:
                    seen.add(x)
                    items.append(x)

    items.sort()
    idx = {name: i for i, name in enumerate(items)}
    n = len(items)
    if n == 0:
        sys.exit("ERROR: no pairs parsed from mash dist output.")
    for a in items:
        for b in items:
            if (a, b) not in pairs and (b, a) not in pairs:
                sys.exit(f"ERROR: missing Mash pair {a}, {b}; cannot impute zero distance.")
            if a == b and pairs.get((a, b)) != 0:
                sys.exit("ERROR: nonzero Mash self-distance.")
    mat = [[0.0] * n for _ in range(n)]
    for (a, b), d in pairs.items():
        i, j = idx[a], idx[b]
        mat[i][j] = d
        mat[j][i] = d  # enforce symmetry

    with open(args.out, "w") as out:
        out.write("item\t" + "\t".join(items) + "\n")
        for i, name in enumerate(items):
            out.write(name + "\t" + "\t".join(f"{mat[i][j]:.6f}" for j in range(n)) + "\n")

    sys.stderr.write(f"[mash_to_matrix] wrote {n}x{n} matrix -> {args.out}\n")


if __name__ == "__main__":
    main()
