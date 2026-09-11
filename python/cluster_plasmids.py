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


def effective_threshold(threshold, len_a, len_b, size_correction_per_percent=0.0, size_correction_cap_pct=40.0):
    """Loosen (never tighten) the threshold as a plasmid pair's size difference grows, per
    Scherff et al. (Microbiol Spectrum 2024, 10.1128/spectrum.02100-24): effective threshold =
    threshold + size_correction_per_percent * min(size_diff_pct, size_correction_cap_pct), where
    size_diff_pct = 100 * |len_a - len_b| / max(len_a, len_b). Disabled (returns threshold
    unchanged) when size_correction_per_percent is 0 or either length is falsy/unknown."""
    if not size_correction_per_percent or not len_a or not len_b:
        return threshold
    size_diff_pct = 100.0 * abs(len_a - len_b) / max(len_a, len_b)
    return threshold + size_correction_per_percent * min(size_diff_pct, size_correction_cap_pct)


def cluster_single(items, dist, threshold, lengths=None, size_correction_per_percent=0.0, size_correction_cap_pct=40.0):
    n = len(items)
    uf = UF(n)
    for i in range(n):
        for j in range(i + 1, n):
            t = threshold if lengths is None else effective_threshold(
                threshold, lengths[i], lengths[j], size_correction_per_percent, size_correction_cap_pct)
            if dist[i][j] <= t:
                uf.union(i, j)
    comp = {}
    for i in range(n):
        comp.setdefault(uf.find(i), []).append(i)
    return list(comp.values())


def cluster_complete(items, dist, threshold, lengths=None, size_correction_per_percent=0.0, size_correction_cap_pct=40.0):
    """Greedy complete-linkage: agglomerate only if the pairwise max stays <= t."""
    n = len(items)
    clusters = [[i] for i in range(n)]

    def can_merge(a, b):
        if lengths is None:
            return all(dist[x][y] <= threshold for x in a for y in b)
        return all(dist[x][y] <= effective_threshold(
            threshold, lengths[x], lengths[y], size_correction_per_percent, size_correction_cap_pct)
            for x in a for y in b)

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


def read_lengths(path, items):
    """TSV with plasmid_id and length columns (e.g. plasmid_index.tsv); header-based lookup
    since such files may carry many other columns in no fixed order."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "plasmid_id" not in header or "length" not in header:
            sys.exit("ERROR: --lengths file must have plasmid_id and length columns.")
        pid_i, len_i = header.index("plasmid_id"), header.index("length")
        lengths = {}
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            lengths[f[pid_i]] = f[len_i]
    missing = sorted(it for it in items if it not in lengths)
    if missing:
        sys.exit(f"ERROR: --lengths file is missing entries for: {missing[:5]}{'...' if len(missing) > 5 else ''}")
    resolved = {}
    for it in items:
        try:
            value = float(lengths[it])
        except ValueError:
            sys.exit(f"ERROR: --lengths value for '{it}' is not a number: {lengths[it]!r}")
        if not math.isfinite(value) or value <= 0:
            # A silently-ignored bad length would make effective_threshold() fall back to the
            # uncorrected threshold for only this item's pairs, with no signal anywhere that
            # size correction was partially disabled -- fail loudly instead.
            sys.exit(f"ERROR: --lengths value for '{it}' must be a positive number, got {lengths[it]!r}")
        resolved[it] = value
    return [resolved[it] for it in items]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--linkage", choices=["single", "complete"], default="single")
    ap.add_argument("--out", required=True)
    ap.add_argument("--lengths", help="TSV with plasmid_id and length columns (e.g. plasmid_index.tsv); required to enable size correction")
    ap.add_argument("--size-correction-per-percent", type=float, default=0.0,
                     help="Loosen the threshold by this amount per 1%% pairwise length difference (0 = disabled); see Scherff et al. 2024, doi:10.1128/spectrum.02100-24")
    ap.add_argument("--size-correction-cap-pct", type=float, default=40.0,
                     help="Size-difference percentage beyond which the correction stops growing")
    args = ap.parse_args()
    if not math.isfinite(args.threshold) or not 0 <= args.threshold <= 1:
        ap.error("threshold must be finite and within [0, 1]")
    if args.size_correction_per_percent < 0:
        ap.error("size-correction-per-percent must be >= 0")
    if not 0 < args.size_correction_cap_pct <= 100:
        ap.error("size-correction-cap-pct must be in (0, 100]")
    if args.size_correction_per_percent and not args.lengths:
        ap.error("--lengths is required when --size-correction-per-percent is nonzero")

    items, idx, dist = read_matrix(args.matrix)
    lengths = read_lengths(args.lengths, items) if args.lengths else None
    if args.linkage == "single":
        comps = cluster_single(items, dist, args.threshold, lengths, args.size_correction_per_percent, args.size_correction_cap_pct)
    else:
        comps = cluster_complete(items, dist, args.threshold, lengths, args.size_correction_per_percent, args.size_correction_cap_pct)

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
    # under single linkage, exposing chained membership beyond the direct threshold. Always
    # measured against the base --threshold, even when size correction is enabled: a single
    # per-unit margin cannot represent per-pair effective thresholds that vary by member length.
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
    correction_note = f", size-corrected +{args.size_correction_per_percent}/1% up to {args.size_correction_cap_pct}%" if args.size_correction_per_percent else ""
    sys.stderr.write(
        f"[cluster] {len(items)} plasmids -> {n_clusters} plasmid units "
        f"({singletons} singletons) at threshold {args.threshold} "
        f"({args.linkage} linkage{correction_note})\n"
    )


if __name__ == "__main__":
    main()
