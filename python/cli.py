"""Portable, validated core runner for reconstructed plasmid inputs."""
import argparse
import json
import math
import platform
import random
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from . import __version__, PROJECT_NAME, PROJECT_TITLE, PROJECT_TAGLINE
from .contracts import (MODES, INDEX_FIELDS, FEATURE_FIELDS, checksum, metadata,
                        normalize_features, read_tsv, validate_manifest, write_tsv)
from .cluster_plasmids import cluster_single, cluster_complete, read_matrix
from .report import build_report
from .quality import assess, load_evidence, QUALITY_FIELDS
from .external_typing import EXTERNAL_TYPING_FIELDS, load_external_typing
from .containment import CONTAINMENT_FIELDS, detect_containment
from .multilayer_network import MULTILAYER_EDGE_FIELDS, multilayer_edges, write_multilayer_graphml


def step(script, *args):
    subprocess.run([sys.executable, str(Path(__file__).with_name(script)), *map(str, args)], check=True)


def aggregate(index, features, meta, assignments, statuses=()):
    members, groups = defaultdict(list), defaultdict(list)
    by_meta = {r["isolate_id"]: r for r in meta}
    for row in index:
        if row["plasmid_id"] in assignments:
            members[assignments[row["plasmid_id"]]].append(row)
    for row in features:
        # Category + canonical gene/product label is descriptive, NOT inferred orthology.
        key = (row["plasmid_unit"], row.get("functional_category", ""),
               row.get("amr_gene") or row.get("gene_symbol") or row.get("product_name") or row["feature_id"])
        groups[key].append(row)
    output = []
    evaluated = {(r['plasmid_id'], r['stage']) for r in statuses if r['status'] == 'complete'}
    for (unit, category, label), rows in sorted(groups.items()):
        carrier_plasmids = {r["plasmid_id"] for r in rows}
        carrier_isolates = {r["isolate_id"] for r in rows}
        meta_rows = [by_meta[i] for i in carrier_isolates]
        dates = sorted(r["date"] for r in meta_rows if r.get("date"))
        prevalence = len(carrier_plasmids) / len(members[unit])
        # Completion follows the actual caller, not the assigned product category.
        stages = {'amr' if r.get('annotation_engine') == 'amrfinder' else 'genes'
                  for r in rows if r.get('annotation_engine') in {'amrfinder', 'bakta', 'prokka'}}
        n_evaluated = sum(bool(stages) and all((r['plasmid_id'], stage) in evaluated for stage in stages)
                          for r in members[unit])
        output.append(dict(plasmid_unit=unit, feature_id=label, functional_category=category,
                           prevalence_in_unit=prevalence,
                           core_or_accessory="observed_in_all" if prevalence == 1 else "observed_in_subset",
                           n_plasmids=len(carrier_plasmids), denominator_plasmids=len(members[unit]),
                           n_evaluated_plasmids=n_evaluated, annotation_coverage=n_evaluated/len(members[unit]),
                           n_isolates=len(carrier_isolates), n_organisms=len({r['organism'] for r in meta_rows if r.get('organism')}),
                           n_locations=len({r['location'] for r in meta_rows if r.get('location')}),
                           first_date=dates[0] if dates else "", last_date=dates[-1] if dates else "",
                           representative_sequence=min(carrier_plasmids),
                           completeness="unresolved", interpretation="annotation-label prevalence; no validated orthology/core-gene inference"))
    return output


def run(args):
    if not math.isfinite(args.threshold) or not 0 <= args.threshold <= 1:
        raise ValueError("threshold must be finite and in [0,1]")
    if args.k < 1 or args.min_length < args.k or args.sketch_size < 1:
        raise ValueError("Require k >= 1, min-length >= k, and sketch-size >= 1")
    containment_min_ratio = getattr(args, 'containment_min_ratio', 0.5)
    containment_max_ratio = getattr(args, 'containment_max_ratio', 0.95)
    containment_max_distance = getattr(args, 'containment_max_distance', None)
    if containment_max_distance is None:
        containment_max_distance = args.threshold
    if not 0 < containment_min_ratio < containment_max_ratio <= 1:
        raise ValueError("Require 0 < containment-min-ratio < containment-max-ratio <= 1")
    if not 0 <= containment_max_distance <= 1:
        raise ValueError("containment-max-distance must be in [0,1]")
    meta = metadata(args.metadata)
    index, sequences = validate_manifest(args.manifest, meta, args.mode, args.min_length)
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory is not empty: {out}. Use a fresh directory to avoid stale results.")
    # Validate feature records before doing any expensive work.
    normalize_features(args.features, index, sequences, {})
    if args.features and getattr(args, 'annotation_config', None):
        raise ValueError("Use either imported --features or --annotation-config to avoid double-counting annotation sources")
    out.mkdir(parents=True, exist_ok=True)
    record = dict(schema_version="1.0", framework_version=__version__, status="running",
                  project_name=PROJECT_NAME, project_title=PROJECT_TITLE, project_tagline=PROJECT_TAGLINE,
                  started_at=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
                  mode=args.mode, dataset_kind="synthetic demonstration" if getattr(args, "synthetic", False) else "user supplied",
                  distance_engine=args.engine, threshold=args.threshold, linkage=args.linkage,
                  k=args.k, sketch_size=args.sketch_size if args.engine == "mash" else None,
                  min_length=args.min_length,
                  manifest_sha256=checksum(args.manifest), metadata_sha256=checksum(args.metadata),
                  features_sha256=checksum(args.features) if args.features else None,
                  annotation_status="imported; absence not evaluated" if args.features else "not_evaluated",
                  quality_policy="evidence-based research rules v1.0; tiers are not calibrated probabilities",
                  sequence_checksum_definition="sha256 of uppercase contig sequences joined by newline, no terminal newline; input contig order",
                  plasbench_dependency=False)
    provenance = out / "run_provenance.json"
    provenance.write_text(json.dumps(record, indent=2), encoding="utf-8")
    try:
        annotation_statuses, typing, features_path = [], {}, args.features
        if getattr(args, 'annotation_config', None):
            from .annotation import annotate_candidates
            features_path, annotation_statuses, typing = annotate_candidates(index, sequences, args.annotation_config,
                out, getattr(args, 'annotation_cache', None), getattr(args, 'threads', 1), getattr(args, 'annotation_profile', 'essential'))
            record['annotation_status'] = 'automatic; per-candidate stage completion recorded'
            record['annotation_config_sha256'] = checksum(args.annotation_config)
        evidence = load_evidence(getattr(args, 'quality_evidence', None), index)
        quality_rows = assess(index, sequences, typing, evidence)
        write_tsv(out/'biological_quality.tsv', QUALITY_FIELDS, quality_rows)
        if getattr(args, 'quality_evidence', None):
            record['quality_evidence_sha256'] = checksum(args.quality_evidence)
        # External typing/taxonomy cross-references are opaque identifiers, reported but never
        # passed to assess() as `typing` -- they must not influence quality/confidence tiers.
        external_typing = load_external_typing(getattr(args, 'external_typing', None), index)
        if getattr(args, 'external_typing', None):
            record['external_typing_sha256'] = checksum(args.external_typing)
        accepted = [r for r in index if r["quality_status"] != "rejected"]
        if external_typing:
            write_tsv(out / "typing_crossreference.tsv", EXTERNAL_TYPING_FIELDS,
                      [external_typing[r["plasmid_id"]] for r in accepted if r["plasmid_id"] in external_typing])
        write_tsv(out / "validation.tsv", INDEX_FIELDS, index)
        pdir = out / "plasmids"
        pdir.mkdir()
        for row in accepted:
            file = pdir / (row["plasmid_id"] + ".fasta")
            file.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences[row["plasmid_id"]]), encoding="utf-8")
        # Output index uses portable relative paths; validation retains original input provenance.
        write_tsv(out / "plasmid_index.tsv", INDEX_FIELDS,
                  [dict(r, fasta_path=f"plasmids/{r['plasmid_id']}.fasta") for r in accepted])
        fields = list(dict.fromkeys(["isolate_id", "chromosomal_cluster", "date", "location", "organism"] + [k for r in meta for k in r]))
        write_tsv(out / "metadata.tsv", fields, meta)
        matrix = out / "plasmid_matrix.tsv"
        if not accepted:
            matrix.write_text("item\n", encoding="utf-8")
            record["distance_status"] = "no_accepted_plasmids"
        elif args.engine == "kmer":
            if len(accepted) > 200:
                raise ValueError("The exact k-mer engine is limited to 200 plasmids; use Mash for larger cohorts")
            step("kmer_distance.py", "--dir", pdir, "--k", args.k, "--out", matrix)
            record["distance_status"] = "computed"
        else:
            mash = shutil.which("mash")
            if not mash:
                raise ValueError("Mash is not on PATH. Install env/environment.yml or explicitly use --engine kmer for small sets")
            record["mash_version"] = subprocess.check_output([mash, "--version"], text=True).strip()
            commands = [[mash, "sketch", "-k", str(args.k), "-s", str(args.sketch_size), "-o", str(out / "plasmids_sketch"),
                         *[str(pdir / (r["plasmid_id"] + ".fasta")) for r in accepted]],
                        [mash, "dist", str(out / "plasmids_sketch.msh"), str(out / "plasmids_sketch.msh")]]
            record["commands"] = commands
            with (out / "mash.log").open("w", encoding="utf-8") as log:
                subprocess.run(commands[0], stdout=log, stderr=log, check=True)
                with (out / "mash_dist.tsv").open("w", encoding="utf-8") as dest:
                    subprocess.run(commands[1], stdout=dest, stderr=log, check=True)
            step("mash_to_matrix.py", "--dist", out / "mash_dist.tsv", "--out", matrix)
            record["distance_status"] = "computed"
        step("cluster_plasmids.py", "--matrix", matrix, "--threshold", args.threshold, "--linkage", args.linkage, "--out", out / "plasmid_clusters.tsv")
        assignments = {r["plasmid_id"]: r["plasmid_unit"] for r in read_tsv(out / "plasmid_clusters.tsv")}
        step("build_network.py", "--index", out / "plasmid_index.tsv", "--clusters", out / "plasmid_clusters.tsv",
             "--metadata", out / "metadata.tsv", "--out-prefix", out / "network")
        step("discordance.py", "--isolate-units", out / "network.isolate_units.tsv", "--metadata", out / "metadata.tsv", "--out-prefix", out / "discordance")
        features = normalize_features(features_path, index, sequences, assignments)
        write_tsv(out / "functional_features.tsv", FEATURE_FIELDS, features)
        functions = aggregate(accepted, features, meta, assignments, annotation_statuses)
        fields = "plasmid_unit feature_id functional_category prevalence_in_unit core_or_accessory n_plasmids denominator_plasmids n_evaluated_plasmids annotation_coverage n_isolates n_organisms n_locations first_date last_date representative_sequence completeness interpretation".split()
        write_tsv(out / "plasmid_unit_function.tsv", fields, functions)
        ids, _, distances = read_matrix(matrix)
        if set(ids) != {r["plasmid_id"] for r in accepted}:
            raise ValueError("Matrix IDs do not match imported plasmids")
        method = cluster_single if args.linkage == "single" else cluster_complete
        sweep = []
        for threshold in sorted({0.0, args.threshold / 2, args.threshold, min(1.0, args.threshold * 2)}):
            comps = method(ids, distances, threshold)
            sweep.append(dict(threshold=threshold, linkage=args.linkage, n_units=len(comps), n_singletons=sum(len(c) == 1 for c in comps)))
        write_tsv(out / "threshold_sensitivity.tsv", ["threshold", "linkage", "n_units", "n_singletons"], sweep)
        # Distinguish direct threshold support from links induced by single-linkage chains.
        pos = {pid: i for i, pid in enumerate(ids)}
        members = defaultdict(list)
        for row in accepted:
            members[(row["isolate_id"], assignments[row["plasmid_id"]])].append(row["plasmid_id"])
        edge_details = []
        eligible = defaultdict(set)
        for f in features:
            if f["headline_eligible"] == "true":
                eligible[f["plasmid_id"]].add(f["amr_gene"])
        for edge in read_tsv(out / "network.edges.tsv"):
            for unit in edge["shared_units"].split(","):
                left, right = members[(edge["source"], unit)], members[(edge["target"], unit)]
                distance = min(distances[pos[a]][pos[b]] for a in left for b in right)
                shared_args = set().union(*(eligible[p] for p in left)) & set().union(*(eligible[p] for p in right))
                edge_details.append(dict(source=edge["source"], target=edge["target"], plasmid_unit=unit,
                                         minimum_distance=distance, direct_threshold_support=str(distance <= args.threshold).lower(),
                                         threshold_margin=args.threshold - distance,
                                         shared_args=";".join(sorted(shared_args)),
                                         interpretation="candidate sharing link; direct transmission unproven"))
        write_tsv(out / "network.edge_evidence.tsv", "source target plasmid_unit minimum_distance direct_threshold_support threshold_margin shared_args interpretation".split(), edge_details)
        if getattr(args, 'multilayer_network', False):
            meta_by_iso = {r["isolate_id"]: r for r in meta}
            layered = multilayer_edges(edge_details, meta_by_iso)
            write_tsv(out / "network.multilayer_edges.tsv", MULTILAYER_EDGE_FIELDS, layered)
            write_multilayer_graphml(out / "network.multilayer.graphml", layered, meta_by_iso.keys())
        containment_rows = detect_containment(accepted, ids, pos, distances,
            containment_min_ratio, containment_max_ratio, containment_max_distance)
        write_tsv(out / "containment_candidates.tsv", CONTAINMENT_FIELDS, containment_rows)
        record.update(status="complete", completed_at=datetime.now(timezone.utc).isoformat(),
                      n_isolates=len(meta), n_plasmids=len(accepted), n_rejected=len(index) - len(accepted), n_units=len(set(assignments.values())))
        # Report generation must succeed before the run is recorded as complete on disk.
        build_report(out, record, index, meta, assignments, features, functions, edge_details, sequences)
        print(f"Complete: {out / 'REPORT.html'}")
    except BaseException as error:
        record.update(status="failed", error=str(error))
        provenance.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        try:
            from .report_output import failure_report
            failure_report(out, record, index, meta, sequences)
        except Exception as report_error:
            record["report_error"] = str(report_error)
            provenance.write_text(json.dumps(record, indent=2), encoding="utf-8")
        raise


def demo(args):
    base = args.out.resolve()
    if base.exists() and any(base.iterdir()):
        raise ValueError("Demo output must be a fresh directory")
    inputs = base / "inputs"
    inputs.mkdir(parents=True)
    rng = random.Random(42)
    a = "".join(rng.choice("ACGT") for _ in range(1200))
    b = "".join(rng.choice("ACGT") for _ in range(900))
    def mutate(seq):
        return "".join(rng.choice([c for c in "ACGT" if c != base]) if rng.random() < 0.01 else base for base in seq)
    rows, features = [], []
    for iso, name, seq in [("iso1", "pA", a), ("iso2", "pA", a), ("iso2", "pB", b), ("iso3", "pA", mutate(a)), ("iso4", "pB", mutate(b))]:
        pid = f"{iso}__{name}"
        (inputs / f"{pid}.fasta").write_text(f">{pid}\n{seq}\n", encoding="utf-8")
        rows.append(dict(isolate_id=iso, plasmid_id=pid, fasta_path=f"{pid}.fasta", source_tool="synthetic_fixture", source_tool_version="1"))
        genes = [("demo_ARG", "amr", 40, 210), ("demo_rep", "replication", 250, 420)] if name == "pA" else [("demo_merA", "metal_resistance", 30, 210)]
        for gene, category, start, end in genes:
            features.append(dict(plasmid_id=pid, feature_id=f"{pid}_{gene}", gene_symbol=gene,
                                 product_name="Illustrative annotation; not detected from sequence", start=start, end=end,
                                 strand="+" if category == "amr" else "-", functional_category=category,
                                 amr_gene=gene if category == "amr" else "", drug_class="demo_class" if category == "amr" else "",
                                 hit_class="curated", annotation_confidence="high", annotation_engine="synthetic_fixture",
                                 database_name="synthetic_not_a_biological_database", database_version="1"))
    write_tsv(inputs / "manifest.tsv", ["isolate_id", "plasmid_id", "fasta_path", "source_tool", "source_tool_version"], rows)
    write_tsv(inputs / "features.tsv", FEATURE_FIELDS, features)
    write_tsv(inputs / "metadata.tsv", ["isolate_id", "chromosomal_cluster", "location", "date", "organism"],
              [dict(isolate_id=f"iso{i}", chromosomal_cluster=f"CC{1 if i < 3 else 2 if i < 5 else 3}",
                    location=f"Synthetic site {1 if i < 3 else 2}", date=f"2026-01-{i:02}", organism="Synthetic organism") for i in range(1, 6)])
    run(argparse.Namespace(manifest=inputs / "manifest.tsv", metadata=inputs / "metadata.tsv", features=inputs / "features.tsv",
                           out=base / "results", mode="precomputed", min_length=200, threshold=0.05,
                           linkage="single", k=15, sketch_size=1000, engine=args.engine, synthetic=True))


def main():
    parser = argparse.ArgumentParser(description=PROJECT_TITLE, epilog=PROJECT_TAGLINE)
    parser.add_argument("--version", action="version", version=f"{PROJECT_NAME} {__version__}")
    subs = parser.add_subparsers(dest="command", required=True)
    p = subs.add_parser("run", help="Validate reconstructed plasmids and build an offline evidence report")
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--features", type=Path, help="Normalized functional TSV, 1-based inclusive coordinates")
    p.add_argument("--annotation-config", type=Path, help="Run local gene, AMR and mobility annotation using a versioned JSON config")
    p.add_argument("--annotation-profile", choices=('essential', 'custom'), default='essential')
    p.add_argument("--annotation-cache", type=Path)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--quality-evidence", type=Path, help="Checksum-linked classification/read/contamination evidence TSV")
    p.add_argument("--external-typing", type=Path, help="Checksum-linked cross-reference TSV for opaque external identifiers (e.g. MOB-suite cluster IDs, COPLA PTU); never influences quality tiers")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--mode", choices=MODES, default="precomputed")
    p.add_argument("--engine", choices=("mash", "kmer"), default="mash")
    p.add_argument("--threshold", type=float, required=True, help="Study-selected distance threshold; not a transmission probability")
    p.add_argument("--linkage", choices=("single", "complete"), default="complete")
    p.add_argument("--k", type=int, default=21)
    p.add_argument("--sketch-size", type=int, default=10000)
    p.add_argument("--min-length", type=int, default=200)
    p.add_argument("--containment-min-ratio", type=float, default=0.5, help="Minimum small/large length ratio for the containment heuristic")
    p.add_argument("--containment-max-ratio", type=float, default=0.95, help="Maximum small/large length ratio for the containment heuristic (equal-length pairs are excluded)")
    p.add_argument("--containment-max-distance", type=float, help="Maximum pairwise distance for the containment heuristic (defaults to --threshold)")
    p.add_argument("--multilayer-network", action="store_true", help="Also write a per-plasmid-unit multilayer GraphML/TSV, each edge tagged with cluster_relation")
    p.set_defaults(func=run)
    p = subs.add_parser("demo", help="Run seeded synthetic data with illustrative functional annotations")
    p.add_argument("--out", type=Path, default=Path("results_demo"))
    p.add_argument("--engine", choices=("kmer", "mash"), default="kmer")
    p.set_defaults(func=demo)
    p = subs.add_parser("check", help="Read-only runtime check")
    p.set_defaults(func=lambda args: print(json.dumps({"python": platform.python_version(), "standalone_core": "available", **{t: shutil.which(t) for t in ('mash','mob_recon','bakta','prokka','amrfinder','mob_typer')}, "automatic_annotation": "available; configure versioned local databases with --annotation-config"}, indent=2)))
    from .importers import import_inputs, ADAPTERS
    p = subs.add_parser('import', help='Import native tool outputs into the common validated manifest')
    p.add_argument('adapter', choices=ADAPTERS)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--input', type=Path)
    g.add_argument('--samples', type=Path, help='TSV: isolate_id, input_path, optional record_ids and metadata')
    p.add_argument('--isolate-id')
    p.add_argument('--record-ids', help='Comma-separated explicit plasmid contig selections')
    p.add_argument('--group-contigs', action='store_true', help='Treat all records as ONE candidate (generic/PlasBench only)')
    p.add_argument('--tool-version', default='unreported')
    p.add_argument('--out', type=Path, required=True)
    p.set_defaults(func=import_inputs)
    from .calibration import add_parser as calibration_parser
    calibration_parser(subs)
    from .report_output import add_parser as report_parser
    report_parser(subs)
    from .population_summary import add_parser as population_summary_parser
    population_summary_parser(subs)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        parser.exit(1, f"ERROR: {error}\n")
