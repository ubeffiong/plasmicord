#!/usr/bin/env python3
"""
build_network.py -- build the outbreak-relevant ISOLATE-level plasmid-sharing
network from plasmid-unit assignments.

Two isolates are linked if they carry plasmids assigned to the same plasmid unit
(PU). Edge weight = number of distinct PUs the two isolates share. This is the
network a public health lab reads: connected isolates may be linked by plasmid
transmission even when their chromosomes are unrelated.

Inputs
------
--index    : plasmid_index.tsv  (columns: plasmid_id, isolate_id, ...)
--clusters : plasmid_clusters.tsv (columns: plasmid_id, plasmid_unit)
--metadata : metadata.tsv (must contain isolate_id column; used to include ALL
             isolates as nodes, even plasmid-free ones, and to annotate nodes
             with chromosomal_cluster/date/location if present)
--out-prefix : writes
     <prefix>.isolate_units.tsv   isolate_id -> comma-separated PUs
     <prefix>.edges.tsv           source, target, weight, shared_units
     <prefix>.graphml             network for Cytoscape / Gephi / microreact

Standard library only.
"""

import argparse
import sys
from collections import defaultdict
from xml.sax.saxutils import escape as xml_escape


def escape(value):
    return xml_escape(value, {'"': '&quot;', "'": '&apos;'})


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            # pad short rows
            while len(f) < len(header):
                f.append("")
            rows.append(dict(zip(header, f)))
    return header, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    _, index_rows = read_tsv(args.index)
    _, clust_rows = read_tsv(args.clusters)
    meta_header, meta_rows = read_tsv(args.metadata)

    if "isolate_id" not in meta_header:
        sys.exit("ERROR: metadata.tsv must have an 'isolate_id' column.")
    for rows, field in ((meta_rows, "isolate_id"), (index_rows, "plasmid_id"), (clust_rows, "plasmid_id")):
        values = [r.get(field, "") for r in rows]
        if any(not v for v in values) or len(values) != len(set(values)):
            sys.exit(f"ERROR: {field} must be unique and non-empty in its table.")
    known_isolates = {r["isolate_id"] for r in meta_rows}
    if any(r["isolate_id"] not in known_isolates for r in index_rows):
        sys.exit("ERROR: plasmid index references isolates missing from metadata.")
    if {r["plasmid_id"] for r in index_rows} != {r["plasmid_id"] for r in clust_rows}:
        sys.exit("ERROR: index and clusters must contain exactly the same plasmids.")

    plasmid_to_isolate = {r["plasmid_id"]: r["isolate_id"] for r in index_rows}
    plasmid_to_unit = {r["plasmid_id"]: r["plasmid_unit"] for r in clust_rows}

    # isolate -> set of PUs
    isolate_units = defaultdict(set)
    for pid, iso in plasmid_to_isolate.items():
        pu = plasmid_to_unit.get(pid)
        if pu:
            isolate_units[iso].add(pu)

    # Node set = every isolate in metadata (include plasmid-free ones).
    all_isolates = [r["isolate_id"] for r in meta_rows]
    meta_by_iso = {r["isolate_id"]: r for r in meta_rows}
    for iso in all_isolates:
        isolate_units.setdefault(iso, set())

    # PU -> isolates carrying it
    unit_to_isolates = defaultdict(set)
    for iso, pus in isolate_units.items():
        for pu in pus:
            unit_to_isolates[pu].add(iso)

    # Edges: for each PU, connect all isolate pairs carrying it.
    edge_units = defaultdict(set)  # (a,b) sorted -> set of shared PUs
    for pu, isos in unit_to_isolates.items():
        isos = sorted(isos)
        for i in range(len(isos)):
            for j in range(i + 1, len(isos)):
                edge_units[(isos[i], isos[j])].add(pu)

    # ---- write isolate_units.tsv ----
    with open(args.out_prefix + ".isolate_units.tsv", "w") as fh:
        fh.write("isolate_id\tn_plasmid_units\tplasmid_units\n")
        for iso in all_isolates:
            pus = sorted(isolate_units[iso])
            fh.write(f"{iso}\t{len(pus)}\t{','.join(pus)}\n")

    # ---- write edges.tsv ----
    with open(args.out_prefix + ".edges.tsv", "w") as fh:
        fh.write("source\ttarget\tweight\tshared_units\n")
        for (a, b), pus in sorted(edge_units.items()):
            fh.write(f"{a}\t{b}\t{len(pus)}\t{','.join(sorted(pus))}\n")

    # ---- write GraphML ----
    extra_attrs = [c for c in ("chromosomal_cluster", "date", "location", "organism")
                   if c in meta_header]
    write_graphml(args.out_prefix + ".graphml", all_isolates, meta_by_iso,
                  isolate_units, edge_units, extra_attrs)

    n_edges = len(edge_units)
    n_connected = len({x for e in edge_units for x in e})
    sys.stderr.write(
        f"[network] {len(all_isolates)} isolates, {n_edges} plasmid-sharing edges, "
        f"{n_connected} isolates in >=1 sharing link.\n"
    )


def write_graphml(path, isolates, meta_by_iso, isolate_units, edge_units, extra_attrs):
    # Define keys
    keys = [
        ('d_units', 'node', 'n_plasmid_units', 'int'),
        ('d_unitlist', 'node', 'plasmid_units', 'string'),
    ]
    for i, a in enumerate(extra_attrs):
        keys.append((f'd_meta{i}', 'node', a, 'string'))
    keys.append(('e_weight', 'edge', 'weight', 'int'))
    keys.append(('e_shared', 'edge', 'shared_units', 'string'))
    attr_key = {a: f'd_meta{i}' for i, a in enumerate(extra_attrs)}

    lines = []
    lines.append('<?xml version="1.0" encoding="UTF-8"?>')
    lines.append('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">')
    for kid, dom, name, typ in keys:
        lines.append(
            f'  <key id="{kid}" for="{dom}" attr.name="{escape(name)}" attr.type="{typ}"/>'
        )
    lines.append('  <graph edgedefault="undirected">')

    for iso in isolates:
        pus = sorted(isolate_units.get(iso, set()))
        lines.append(f'    <node id="{escape(iso)}">')
        lines.append(f'      <data key="d_units">{len(pus)}</data>')
        lines.append(f'      <data key="d_unitlist">{escape(",".join(pus))}</data>')
        for a in extra_attrs:
            val = meta_by_iso.get(iso, {}).get(a, "")
            lines.append(f'      <data key="{attr_key[a]}">{escape(val)}</data>')
        lines.append('    </node>')

    eid = 0
    for (a, b), shared in sorted(edge_units.items()):
        lines.append(f'    <edge id="e{eid}" source="{escape(a)}" target="{escape(b)}">')
        lines.append(f'      <data key="e_weight">{len(shared)}</data>')
        lines.append(f'      <data key="e_shared">{escape(",".join(sorted(shared)))}</data>')
        lines.append('    </edge>')
        eid += 1

    lines.append('  </graph>')
    lines.append('</graphml>')
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
