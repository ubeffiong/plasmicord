"""Self-contained, offline HTML evidence report and machine-readable payload."""
from .contracts import read_tsv

SAFEGUARDS = [
    "Shared plasmid units and functional similarity do not establish direct person-to-person transmission or direction.",
    "Mash-style distance is an approximate k-mer distance, not measured alignment identity or a transmission probability.",
    "Single linkage can connect plasmids through intermediate members; inspect direct threshold support.",
    "Unknown chromosomal clusters are excluded from cross-cluster calls. No detected link does not prove no transmission.",
    "Gene detection does not prove phenotypic resistance or expression; partial ARGs may not be functional.",
    "Broad efflux matches require caution. Plasmid-only analysis can miss chromosomal resistance mutations.",
    "A pathway-associated gene does not establish a complete pathway. Module completeness is unresolved unless separately evaluated.",
    "Headline ARGs require a strict/perfect/curated high-confidence call and engine/database provenance. Loose and partial calls are excluded.",
    "Imported annotations do not establish complete annotation coverage; missing features mean not evaluated, not biological absence.",
    "Biological quality tiers are evidence-based research rules, not calibrated probabilities. Marker evidence alone does not confirm plasmid identity or closure.",
    "Database versions and detection thresholds affect results. Plasmid unit identifiers are run-local and may change when the cohort changes.",
]


def build_report(out, provenance, index, meta, assignments, features, functions, evidence, sequences):
    crosslinks = read_tsv(out / "discordance.crosslinks.tsv")
    edges = read_tsv(out / "network.edges.tsv")
    known = sum(bool(r.get("chromosomal_cluster")) for r in meta)
    annotation_status = read_tsv(out/'annotation_status.tsv') if (out/'annotation_status.tsv').is_file() else []
    quality = read_tsv(out/'biological_quality.tsv') if (out/'biological_quality.tsv').is_file() else []
    amr_evaluated = sum(r.get('stage') == 'amr' and r.get('status') == 'complete' for r in annotation_status)
    counts = dict(isolates=len(meta), plasmids=len(assignments), units=len(set(assignments.values())),
                  sharing_pairs=len(edges), cross_cluster_pairs=len({(r['isolate_a'], r['isolate_b']) for r in crosslinks}),
                  cross_cluster_unit_links=len(crosslinks), typed_isolates=known,
                  observed_arg_plasmids=len({r['plasmid_id'] for r in features if r['headline_eligible'] == 'true'}) if features or amr_evaluated else None)
    payload = dict(schema_version="1.0", provenance=provenance, counts=counts, metadata=meta,
                   plasmids=[dict(r, plasmid_unit=assignments.get(r['plasmid_id'], "")) for r in index],
                   features=features, unit_functions=functions, edges=edges, edge_evidence=evidence,
                   crosslinks=crosslinks, safeguards=SAFEGUARDS, biological_quality=quality, annotation_status=annotation_status,
                   contigs={pid: [{"id": name, "length": len(seq)} for name, seq in records] for pid, records in sequences.items()},
                   sensitivity=read_tsv(out / "threshold_sensitivity.tsv"),
                   plasmid_clusters=read_tsv(out / "plasmid_clusters.tsv") if (out / "plasmid_clusters.tsv").is_file() else [])
    from .report_output import render_dashboard
    return render_dashboard(out, payload)
