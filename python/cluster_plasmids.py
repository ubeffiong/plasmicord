#!/usr/bin/env python3
"""
cluster_plasmids.py -- turn an all-vs-all distance matrix into "plasmid units"
(PUs): reproducible threshold-based clusters, analogous to how cgMLST gives
sequence types.

Method
------
Build a graph where every pair of items with distance <= THRESHOLD is connected,
then take clusters as:
  * single   linkage (default): connected components (union-find). Two plasmids
             join a cluster if EITHER is within threshold of ANY member. Standard
             for Mash-based plasmid clustering (MOB-cluster-like).
  * complete linkage: a plasmid joins a cluster only if it is within threshold of
             EVERY current member (stricter, fewer chimeric clusters).

Cluster IDs are deterministic: clusters are sorted by their lexicographically
smallest member, then numbered PU_0001, PU_0002, ...

Input
-----
--matrix : TSV. First row = header: an empty/ignored first cell, then item ids.
           Each subsequent row = item id, then its distances to every column item.
           Must be square and symmetric (diagonal ~0). Distances are floats.
--threshold : max distance to link two items (e.g. 0.05 for Mash ~= 95% ANI).
--linkage   : single (default) | complete
--out       : TSV written as: item_id <TAB> cluster_id

Standard library only.
"""

import argparse
import math
import sys


def read_matrix(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        items = header[1:]
        if len(items) != len(set(items)) or any(not x for x in items):
            sys.exit("ERROR: matrix identifiers must be unique and non-empty.")
        idx = {name: i for i, name in enumerate(items)}
        n = len(items)
        dist = [[0.0] * n for _ in range(n)]
        row_order = []
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rid = f[0]
            row_order.append(rid)
            if rid not in idx:
                sys.exit(f"ERROR: row id '{rid}' not present in header columns.")
            i = idx[rid]
            vals = f[1:]
            if len(vals) != n:
                sys.exit(f"ERROR: row '{rid}' has {len(vals)} values, expected {n}.")
            for j, v in enumerate(vals):
                dist[i][j] = float(v)
    if sorted(row_order) != sorted(items):
        sys.exit("ERROR: matrix rows and columns must contain the same item ids.")
    for i in range(n):
        for j in range(n):
            value = dist[i][j]
            if not math.isfinite(value) or not 0 <= value <= 1:
                sys.exit("ERROR: distances must be finite and within [0, 1].")
            if abs(value - dist[j][i]) > 1e-9 or (i == j and abs(value) > 1e-9):
                sys.exit("ERROR: matrix must be symmetric with a zero diagonal.")
    return items, idx, dist


# ---- union-find for single linkage ----
class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def cluster_single(items, dist, threshold):
    n = len(items)
    uf = UF(n)
    for i in range(n):
        for j in range(i + 1, n):
            if dist[i][j] <= threshold:
                uf.union(i, j)
    comp = {}
    for i in range(n):
        comp.setdefault(uf.find(i), []).append(i)
    return list(comp.values())


def cluster_complete(items, dist, threshold):
    """Greedy complete-linkage: agglomerate only if the pairwise max stays <= t."""
    n = len(items)
    clusters = [[i] for i in range(n)]

    def can_merge(a, b):
        return all(dist[x][y] <= threshold for x in a for y in b)

    merged = True
    while merged:
        merged = False
        best = None
        best_link = None
        for ci in range(len(clusters)):
            for cj in range(ci + 1, len(clusters)):
                if can_merge(clusters[ci], clusters[cj]):
                    # link distance = max pairwise; prefer smallest such
                    link = max(dist[x][y] for x in clusters[ci] for y in clusters[cj])
                    tie = tuple(sorted(tuple(sorted(items[x] for x in group)) for group in (clusters[ci], clusters[cj])))
                    candidate = (link, tie)
                    if best_link is None or candidate < best_link:
                        best_link = candidate
                        best = (ci, cj)
        if best is not None:
            ci, cj = best
            clusters[ci].extend(clusters[cj])
            del clusters[cj]
            merged = True
    return clusters


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--linkage", choices=["single", "complete"], default="single")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if not math.isfinite(args.threshold) or not 0 <= args.threshold <= 1:
        ap.error("threshold must be finite and within [0, 1]")

    items, idx, dist = read_matrix(args.matrix)
    if args.linkage == "single":
        comps = cluster_single(items, dist, args.threshold)
    else:
        comps = cluster_complete(items, dist, args.threshold)

    # Deterministic ordering + IDs.
    named = []
    for members in comps:
        names = sorted(items[i] for i in members)
        named.append((names, members))
    named.sort(key=lambda pair: pair[0][0])

    width = max(4, len(str(len(named))))
    item_to_cluster = {}
    # unit_threshold_margin: threshold minus the cluster's worst (max) internal pairwise
    # distance. Always >=0 under complete linkage (an algorithm invariant); can go negative
    # under single linkage, exposing chained membership beyond the direct threshold.
    cluster_margin = {}
    for k, (names, members) in enumerate(named, start=1):
        cid = f"PU_{k:0{width}d}"
        for nm in names:
            item_to_cluster[nm] = cid
        if len(members) > 1:
            worst = max(dist[a][b] for a in members for b in members if a != b)
            cluster_margin[cid] = (f"{worst:.6f}", f"{args.threshold - worst:.6f}")
        else:
            cluster_margin[cid] = ("", "")

    with open(args.out, "w") as out:
        out.write("plasmid_id\tplasmid_unit\tunit_max_internal_distance\tunit_threshold_margin\n")
        for nm in items:
            cid = item_to_cluster[nm]
            worst, margin = cluster_margin[cid]
            out.write(f"{nm}\t{cid}\t{worst}\t{margin}\n")

    n_clusters = len(named)
    singletons = sum(1 for names in named if len(names) == 1)
    sys.stderr.write(
        f"[cluster] {len(items)} plasmids -> {n_clusters} plasmid units "
        f"({singletons} singletons) at threshold {args.threshold} "
        f"({args.linkage} linkage)\n"
    )


if __name__ == "__main__":
    main()
