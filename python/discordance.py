#!/usr/bin/env python3
"""
discordance.py -- the headline outbreak-investigation output.

Chromosomal typing (cgMLST / SNP clusters) groups isolates by vertical descent.
Plasmids move horizontally. This step flags where the two disagree:

  * CROSS-LINKS: two isolates that share a plasmid unit but sit in DIFFERENT
    chromosomal clusters. These are candidate plasmid-mediated transmission
    events that pure chromosomal typing would miss -- the reason plasmid-aware
    surveillance matters.

  * (context) same chromosomal cluster but NO shared plasmid unit: isolates that
    are clonally related yet carry different plasmids (plasmid loss/replacement).

Inputs
------
--isolate-units : <prefix>.isolate_units.tsv from build_network.py
                  (isolate_id, n_plasmid_units, plasmid_units)
--metadata      : metadata.tsv with columns isolate_id and chromosomal_cluster
                  (isolates lacking a chromosomal_cluster value are skipped for
                   cross-link calls, with a warning)
--out-prefix    : writes
     <prefix>.crosslinks.tsv    plasmid-mediated links across chromosomal clusters
     <prefix>.summary.txt       counts + a plain-language summary

Standard library only.
"""

import argparse
import sys
from collections import defaultdict
from itertools import combinations


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            while len(f) < len(header):
                f.append("")
            rows.append(dict(zip(header, f)))
    return header, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--isolate-units", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    _, unit_rows = read_tsv(args.isolate_units)
    meta_header, meta_rows = read_tsv(args.metadata)
    if "chromosomal_cluster" not in meta_header:
        sys.exit("ERROR: metadata.tsv needs a 'chromosomal_cluster' column for "
                 "discordance analysis. (Add it, or run only build_network.)")

    chrom = {}
    for r in meta_rows:
        c = r.get("chromosomal_cluster", "").strip()
        chrom[r["isolate_id"]] = c if c.upper() not in ("", "UNKNOWN", "NA", "N/A", "NONE") else None

    # isolate -> set of PUs
    iso_units = {}
    for r in unit_rows:
        pus = set(p for p in r.get("plasmid_units", "").split(",") if p)
        iso_units[r["isolate_id"]] = pus

    # PU -> isolates carrying it
    unit_to_isos = defaultdict(list)
    for iso, pus in iso_units.items():
        for pu in pus:
            unit_to_isos[pu].append(iso)

    n_missing_chrom = sum(1 for iso in iso_units if chrom.get(iso) is None)
    if n_missing_chrom:
        sys.stderr.write(
            f"[discordance] NOTE: {n_missing_chrom} isolate(s) lack a "
            f"chromosomal_cluster; their pairs are excluded from cross-cluster calls.\n"
        )

    # ---- cross-links ----
    crosslinks = []              # (pu, isoA, clA, isoB, clB)
    units_crossing = set()
    for pu, isos in unit_to_isos.items():
        isos = sorted(set(isos))
        for a, b in combinations(isos, 2):
            ca = chrom.get(a) or "UNKNOWN"
            cb = chrom.get(b) or "UNKNOWN"
            if ca != "UNKNOWN" and cb != "UNKNOWN" and ca != cb:
                crosslinks.append((pu, a, ca, b, cb))
                units_crossing.add(pu)

    with open(args.out_prefix + ".crosslinks.tsv", "w") as fh:
        fh.write("plasmid_unit\tisolate_a\tchrom_cluster_a\tisolate_b\tchrom_cluster_b\n")
        for row in sorted(crosslinks):
            fh.write("\t".join(row) + "\n")

    # ---- context: clonal but no shared plasmid ----
    # group isolates by chromosomal cluster; within a cluster, pairs sharing no PU
    by_chrom = defaultdict(list)
    for iso, c in chrom.items():
        if c and iso in iso_units:
            by_chrom[c].append(iso)
    clonal_no_share = 0
    for c, isos in by_chrom.items():
        for a, b in combinations(sorted(isos), 2):
            if not (iso_units[a] & iso_units[b]):
                clonal_no_share += 1

    # ---- summary ----
    n_isolates = len(iso_units)
    n_units = len(unit_to_isos)
    n_shared_units = sum(1 for pu, isos in unit_to_isos.items() if len(set(isos)) > 1)
    with open(args.out_prefix + ".summary.txt", "w") as fh:
        fh.write("Plasmid-aware discordance summary\n")
        fh.write("=================================\n\n")
        fh.write(f"Isolates analysed .......................... {n_isolates}\n")
        fh.write(f"Plasmid units (PUs) ........................ {n_units}\n")
        fh.write(f"PUs shared by >1 isolate ................... {n_shared_units}\n")
        fh.write(f"PUs crossing chromosomal clusters .......... {len(units_crossing)}\n")
        fh.write(f"Cross-cluster plasmid-sharing pairs ........ {len(crosslinks)}\n")
        fh.write(f"Unique cross-cluster isolate pairs .......... {len({(r[1], r[3]) for r in crosslinks})}\n")
        fh.write(f"Isolates with unknown chromosome cluster .... {n_missing_chrom}\n")
        fh.write(f"Clonal pairs (same chrom cluster) w/o shared PU  {clonal_no_share}\n\n")
        if crosslinks:
            fh.write("INTERPRETATION: the cross-cluster pairs above are candidate\n")
            fh.write("plasmid-mediated links that chromosomal typing alone would miss.\n")
            fh.write("Prioritise the plasmid units listed in crosslinks.tsv for review.\n")
        else:
            fh.write("INTERPRETATION: no cross-cluster sharing was detected among\n")
            fh.write("isolates with known clusters. Missing typing limits interpretation.\n")

    sys.stderr.write(
        f"[discordance] {len(units_crossing)} PUs cross chromosomal clusters; "
        f"{len(crosslinks)} cross-cluster isolate pairs. "
        f"See {args.out_prefix}.summary.txt\n"
    )


if __name__ == "__main__":
    main()
