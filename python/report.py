"""Self-contained, offline HTML evidence report and machine-readable payload."""
import html
import json
from pathlib import Path
from . import PROJECT_TITLE, PROJECT_TAGLINE
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
                   sensitivity=read_tsv(out / "threshold_sensitivity.tsv"))
    (out / "report_data.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    data = json.dumps(payload).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    files = ["report_data.json", "run_provenance.json", "validation.tsv", "plasmid_index.tsv", "plasmid_matrix.tsv",
             "plasmid_clusters.tsv", "network.graphml", "network.edges.tsv", "network.edge_evidence.tsv",
             "network.isolate_units.tsv", "discordance.crosslinks.tsv", "discordance.summary.txt",
             "functional_features.tsv", "plasmid_unit_function.tsv", "threshold_sensitivity.tsv"]
    files.extend(name for name in ('biological_quality.tsv','annotation_status.tsv','annotation_provenance.json','mobility_typing.json','plasbench_proteins.tsv') if (out/name).is_file())
    links = "".join(f'<li><a download href="{name}">{name}</a></li>' for name in files)
    template = Path(__file__).with_name("report_template.html").read_text(encoding="utf-8")
    (out / "REPORT.html").write_text(template.replace("__REPORT_DATA__", data).replace("__DOWNLOADS__", links), encoding="utf-8")
    lines = [f"# {PROJECT_TITLE}", "", PROJECT_TAGLINE, "", f"Dataset: **{provenance['dataset_kind']}**. PlasmiCord {provenance['framework_version']}.", "",
             f"{counts['isolates']} isolates; {counts['plasmids']} accepted plasmids; {counts['units']} plasmid units; "
             f"{counts['sharing_pairs']} sharing pairs; {counts['cross_cluster_pairs']} unique cross-cluster pairs "
             f"({counts['cross_cluster_unit_links']} pair-unit links).", "",
             f"Chromosomal typing available for {known}/{len(meta)} isolates. Annotation: {provenance['annotation_status']}.", "",
             f"Engine: {provenance['distance_engine']}; k={provenance['k']}; threshold={provenance['threshold']}; linkage={provenance['linkage']}.", "",
             "Open [REPORT.html](REPORT.html) for the interactive network, functional cargo, gene tracks, quality and downloads.", "",
             "## Interpretation limits", "", *[f"- {s}" for s in SAFEGUARDS], "", "## Evidence files", "", *[f"- [{f}]({f})" for f in files], ""]
    (out / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
