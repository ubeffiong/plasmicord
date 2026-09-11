"""Deterministic report interpretations and auditable offline artifact inventory."""
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
from urllib.parse import quote

RULE_VERSION = "1.0"
PRESENTATION = {"REPORT.html", "REPORT.md", "report_data.json", "run_provenance.json", "output_manifest.json", "REPORT_BUNDLE.zip"}
PREVIEW_LIMIT = 12000
PREVIEW_BUDGET = 1024 * 1024


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def result_files(root):
    """Never traverse symlinks/junctions or expose files outside the result directory."""
    root = Path(root).resolve()
    for base, directories, files in os.walk(root, followlinks=False):
        directories[:] = sorted(d for d in directories if not (Path(base)/d).is_symlink()
                                and (Path(base)/d).resolve().is_relative_to(root))
        for name in sorted(files):
            path = Path(base)/name
            if not path.is_symlink() and path.is_file() and path.resolve().is_relative_to(root):
                yield path


def describe(path):
    name = Path(path).name
    if path.startswith("annotation/"):
        if name == "provenance.json":
            return "Annotation", "Stage command, tool/database identity and cache evidence"
        if name == "tool.log":
            return "Annotation", "Native tool execution log; inspect warnings and failures"
        return "Annotation", "Native annotation intermediate or output; candidate and stage are in its path"
    if path.startswith("calibration/"):
        return "Calibration", "Study-specific training and holdout calibration evidence"
    if path.startswith("plasmids/"):
        return "Candidates", "Normalized accepted candidate FASTA, preserving contig boundaries"
    descriptions = {
        "validation.tsv": ("Quality", "All candidate validation decisions, including exclusions"),
        "biological_quality.tsv": ("Quality", "Sequence, marker, graph and supplied biological evidence"),
        "plasmid_index.tsv": ("Candidates", "Accepted candidate IDs, sequence digests and source metadata"),
        "metadata.tsv": ("Metadata", "Isolate metadata used by this run"),
        "plasmid_matrix.tsv": ("Distances", "Complete pairwise distance matrix"),
        "mash_dist.tsv": ("Distances", "Raw Mash pairwise distance output"),
        "mash.log": ("Distances", "Mash execution log"),
        "plasmid_clusters.tsv": ("Units", "Run-local plasmid-unit membership"),
        "threshold_sensitivity.tsv": ("Units", "Unit/singleton counts across alternative thresholds"),
        "containment_candidates.tsv": ("Units", "Length/similarity heuristic pairs; not alignment-confirmed containment"),
        "typing_crossreference.tsv": ("Annotation", "Opaque external identifiers (e.g. MOB-suite cluster IDs, COPLA PTU); never used to compute quality tiers"),
        "network.multilayer_edges.tsv": ("Sharing", "Per-(isolate pair, plasmid unit) sharing edges tagged with cluster relation; present only with --multilayer-network"),
        "network.multilayer.graphml": ("Sharing", "Multilayer sharing network for Cytoscape or Gephi; present only with --multilayer-network"),
        "population_summary.pu_level.tsv": ("Report", "Cohort-level plasmid-unit aggregation from plasmicord population-summary"),
        "population_summary.metadata_dimension.tsv": ("Report", "Cohort-level metadata-dimension aggregation from plasmicord population-summary"),
        "functional_features.tsv": ("Functions", "Validated contig-coordinate functional evidence"),
        "automatic_features.tsv": ("Annotation", "Normalized automatic-caller output before unit assignment"),
        "plasmid_unit_function.tsv": ("Functions", "Observed function prevalence and annotation coverage"),
        "annotation_status.tsv": ("Annotation", "Candidate-by-stage completion and cache records"),
        "annotation_provenance.json": ("Annotation", "Engine/database versions, hashes and original commands"),
        "mobility_typing.json": ("Annotation", "Candidate-level replicon and mobility predictions"),
        "plasbench_proteins.tsv": ("Functions", "Independent protein export using zero-based half-open coordinates"),
        "run_provenance.json": ("Provenance", "Run state, parameters, input/output hashes and tool versions"),
        "sample_summary.tsv": ("Report", "Per-isolate evidence coverage and descriptive findings"),
        "module_summary.tsv": ("Report", "Module states, observed evidence and next review steps"),
        "interpretations.json": ("Report", "Versioned deterministic findings with rules and evidence links"),
        "report_charts.json": ("Report", "Exact values and denominators behind dashboard charts"),
        "REPORT.html": ("Report", "Interactive offline dashboard"),
        "REPORT.md": ("Report", "Detailed plain-text companion report"),
        "report_data.json": ("Report", "Embedded dashboard data model"),
        "output_manifest.json": ("Provenance", "Final recursive inventory and SHA-256 checksums"),
        "REPORT_BUNDLE.zip": ("Report", "Portable report and all inventoried result files"),
    }
    if name in descriptions:
        return descriptions[name]
    if name.startswith("network."):
        return "Sharing", "Isolate sharing network or pair-unit supporting evidence"
    if name.startswith("discordance."):
        return "Chromosome context", "Known chromosome-cluster discordance evidence"
    return "Other outputs", "Additional file present in the result directory"


def catalog(root, include_generated=False, previews=True):
    root = Path(root).resolve()
    entries, remaining = [], PREVIEW_BUDGET
    text_types = {".tsv", ".csv", ".txt", ".log", ".json", ".jsonl", ".gff", ".gff3", ".fasta", ".fa", ".fna", ".faa", ".ffn", ".md", ".html", ".graphml", ".gfa"}
    for path in result_files(root):
        rel = path.relative_to(root).as_posix()
        if rel in PRESENTATION and not include_generated:
            continue
        if rel in {"REPORT_BUNDLE.zip", "output_manifest.json"}:
            continue
        stat = path.stat()
        module, description = describe(rel)
        entry = dict(path=rel, href="./"+quote(rel, safe="/"), name=path.name, module=module,
                     description=description, size_bytes=stat.st_size, sha256=sha256(path),
                     modified_utc=datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                     media_type=path.suffix.lower().lstrip(".") or "file", preview=None)
        if previews and path.suffix.lower() in text_types and remaining > 0:
            limit = min(PREVIEW_LIMIT, remaining)
            with path.open("rb") as handle:
                data = handle.read(limit)
            if b"\0" not in data:
                entry["preview"] = data.decode("utf-8", errors="replace")
                entry["preview_truncated"] = stat.st_size > len(data)
                remaining -= len(data)
        entries.append(entry)
    return entries


def summary(payload):
    """Findings are rules over measured fields, never inferred transmission probabilities."""
    candidates = payload["plasmids"]
    accepted = [p for p in candidates if p.get("plasmid_unit")]
    metadata = payload["metadata"]
    features = payload["features"]
    statuses = payload.get("annotation_status", [])
    run = payload["provenance"]
    n, total = len(accepted), len(candidates)
    typed = sum(bool(m.get("chromosomal_cluster")) for m in metadata)
    eligible = [f for f in features if f.get("headline_eligible") == "true"]
    observed = {f["plasmid_id"] for f in eligible}
    completed = {(s["plasmid_id"], s["stage"]) for s in statuses if s.get("status") == "complete"}
    accepted_ids = {p["plasmid_id"] for p in accepted}
    completion = {stage: sum((p["plasmid_id"], stage) in completed for p in accepted) for stage in ("genes", "amr", "mobility")}
    qc = Counter(p.get("quality_status") or "uncertain" for p in candidates)
    findings = []

    def finding(rule, state, title, text, evidence, action, section):
        findings.append(dict(rule_id=rule, rule_version=RULE_VERSION, state=state, title=title,
                             interpretation=text, evidence=evidence, next_step=action, section=section))

    failed = run.get("status") != "complete"
    finding("run_state", "fail" if failed else "pass", "Execution state",
            ("The run did not complete. This dashboard describes available intermediate evidence; missing modules are not negative results."
             if failed else "The analysis completed. Completion confirms execution, not biological validity."),
            ["run_provenance.json"], "Inspect the run record and tool logs." if failed else "Review evidence coverage before interpreting links.", "run")
    rejected = qc.get("rejected", 0)
    unresolved = sum(qc.get(k, 0) for k in ("uncertain", "low_confidence"))
    finding("candidate_gate", "warn" if rejected or unresolved else "pass" if total else "not_evaluated",
            "Candidate quality",
            f"{n} candidates have unit assignments out of {total} submitted; {rejected} rejected, {unresolved} low-confidence or uncertain. Quality tiers are heuristic, not validated probabilities.",
            ["validation.tsv", "biological_quality.tsv"], "Inspect exclusion reasons and unresolved sequence/read/closure evidence.", "quality")
    if failed:
        text = "The run is incomplete. No completed chromosome-discordance assessment is claimed; inspect the available intermediate evidence."
    elif not typed:
        text = f"No comparable chromosome clusters were supplied for {len(metadata)} isolates. Chromosomal concordance and discordance cannot be evaluated."
    else:
        text = f"Chromosome typing covers {typed}/{len(metadata)} isolates. {payload['counts']['cross_cluster_pairs']} unique sharing pairs connect distinct known clusters; untyped pairs are excluded."
    finding("typing_coverage", "not_evaluated" if failed or not typed else "warn" if typed < len(metadata) else "pass",
            "Chromosome context", text, ["metadata.tsv", "discordance.crosslinks.tsv"],
            "Supply comparable, namespaced chromosome typing and investigate candidate cross-cluster links.", "sharing")
    indirect = sum(str(e.get("direct_threshold_support")).lower() == "false" for e in payload["edge_evidence"])
    finding("sharing_support", "not_evaluated" if failed else "warn" if indirect else "info", "Sharing and direct support",
            ("No completed sharing summary is available. Intermediate network files, if present, require review with the run failure."
             if failed else f"{len(payload['edges'])} isolate-sharing pairs are represented by {len(payload['edge_evidence'])} pair-unit observations. {indirect} observations lack a direct pairwise distance below the chosen threshold. Neither membership nor direct distance establishes transmission."),
            ["network.edge_evidence.tsv", "plasmid_clusters.tsv"], "Prioritize sequence/structural confirmation and independent epidemiological review.", "sharing")
    clusters = payload.get("plasmid_clusters", [])
    multi_member_units = {c["plasmid_unit"] for c in clusters if c.get("unit_threshold_margin") not in (None, "")}
    chained_units = {c["plasmid_unit"] for c in clusters if c.get("unit_threshold_margin") not in (None, "") and float(c["unit_threshold_margin"]) < 0}
    finding("chaining_risk", "not_evaluated" if failed or not clusters else "warn" if chained_units else "info", "Single-linkage chaining risk",
            ("No completed clustering summary is available; unit threshold-margin diagnostics require the run to finish."
             if failed or not clusters else f"{len(chained_units)}/{len(multi_member_units)} multi-member plasmid units have a negative unit threshold margin: at least one internal member pair exceeds the chosen threshold, so the unit is held together only by chained membership, not mutual similarity within the threshold."),
            ["plasmid_clusters.tsv"], "Inspect unit_max_internal_distance/unit_threshold_margin directly for units flagged this way.", "sharing")
    amr_n = completion["amr"]
    zero_hits = len({pid for pid, stage in completed if stage == "amr" and pid in accepted_ids} - observed)
    finding("amr_coverage", "pass" if n and amr_n == n else "warn" if features or amr_n else "not_evaluated",
            "AMR observations and evaluation coverage",
            f"{len(observed)} accepted plasmids carry {len({f.get('amr_gene') for f in eligible})} distinct eligible ARG labels. AMR caller completion is recorded for {amr_n}/{n} accepted plasmids; {zero_hits} completed candidates have no headline-eligible ARG call. Imported rows alone do not establish search completeness or absence.",
            ["functional_features.tsv", "annotation_status.tsv"], "Inspect caller method, database release, coverage and partial/plus-scope matches before interpreting resistance.", "cargo")
    for stage in ("genes", "mobility"):
        finding(stage+"_coverage", "pass" if n and completion[stage] == n else "warn" if completion[stage] else "not_evaluated",
                stage.capitalize()+" caller coverage", f"{completion[stage]}/{n} accepted plasmids have a completed {stage} stage. Unrun stages remain not evaluated.",
                ["annotation_status.tsv", "annotation_provenance.json"], "Use the configured annotation profile or inspect the recorded external evidence.", "quality")
    missing = {field: sum(not m.get(field) for m in metadata) for field in ("date", "location", "organism")}
    finding("metadata_coverage", "warn" if any(missing.values()) else "pass", "Metadata completeness",
            f"Missing observations: date {missing['date']}, location {missing['location']}, organism {missing['organism']} of {len(metadata)} isolates. Collection timing is descriptive and does not orient transmission.",
            ["metadata.tsv"], "Resolve missing context upstream; do not infer it from plasmid similarity.", "samples")
    sensitivity = payload.get("sensitivity", [])
    counts = [int(s["n_units"]) for s in sensitivity]
    finding("threshold_sensitivity", "warn" if len(set(counts)) > 1 else "info" if counts else "not_evaluated",
            "Threshold sensitivity",
            (f"Unit counts range from {min(counts)} to {max(counts)} across {len(counts)} tested cutoffs. The chosen distance is {run.get('threshold')}; stable counts would not validate the cutoff."
             if counts else "No threshold sweep is available."),
            ["threshold_sensitivity.tsv"], "Compare sensitivity and independent, study-specific calibration.", "distances")
    calibration = payload.get("calibration")
    if calibration:
        c = calibration["record"]
        selected = c.get("selected_threshold")
        finding("calibration", "info" if selected is not None else "not_evaluated", "Attached calibration",
                f"Calibration status: {c.get('status')}; target: {c.get('target')}; training-selected cutoff: {selected}. Current analysis cutoff: {run.get('threshold')}. A report attachment does not change the analysis. {c.get('limitation', '')}",
                ["calibration/calibration.json", "calibration/calibration_metrics.tsv"], "Review training/holdout support and population scope before applying a threshold.", "calibration")
    else:
        finding("calibration", "not_evaluated", "Calibration not attached",
                "No matching calibration result is attached to this report. Threshold sensitivity alone is not calibration.",
                ["threshold_sensitivity.tsv"], "Run plasmicord calibrate and attach its matching result with plasmicord report --calibration.", "calibration")
    duplicates = sum(bool(p.get("duplicate_of")) for p in candidates)
    if duplicates:
        finding("duplicate_sequences", "warn", "Exact duplicate sequences",
                f"{duplicates} candidates duplicate an earlier sequence digest. Biological sharing and technical duplication require separate review.",
                ["validation.tsv"], "Check sample provenance and independent assembly/read evidence.", "quality")
    cached = sum(s.get("cache") == "reused" for s in statuses)
    if cached:
        finding("cache_reuse", "info", "Annotation reuse",
                f"{cached} stage results were reused from checksum-verified cache. Raw files may exist only in the original run; the current stage provenance records the original command and raw-output hash.",
                ["annotation_provenance.json"], "Use the current stage records to distinguish cached results from newly executed calls.", "quality")
    containment = payload.get("containment", [])
    if containment:
        interpretations = Counter(c.get("interpretation") for c in containment)
        finding("containment_candidates", "warn", "Length/similarity containment candidates",
                f"{len(containment)} candidate pairs matched the configured length-ratio and distance heuristic ({dict(interpretations)}). This is not alignment-confirmed containment or a co-integrate/subclone call.",
                ["containment_candidates.tsv"], "Confirm candidate pairs with sequence alignment before treating either as contained within the other.", "sharing")
    typing_crossreference = payload.get("typing_crossreference", [])
    if typing_crossreference:
        finding("typing_coverage_external", "info", "External typing cross-reference",
                f"{len(typing_crossreference)}/{n} accepted plasmids carry an attached external identifier (e.g. MOB-suite cluster, COPLA PTU, PlasmidFinder Inc-type, PLSDB nearest match). These opaque cross-references never influence PlasmiCord's own quality tiers.",
                ["typing_crossreference.tsv"], "Cross-check external calls against their own source database version and confidence before combining with PlasmiCord evidence.", "quality")
    multilayer_edges = payload.get("multilayer_edges", [])
    if multilayer_edges:
        cross = sum(e.get("cluster_relation") == "cross_cluster" for e in multilayer_edges)
        finding("multilayer_network", "info", "Multilayer (isolate pair, plasmid unit) network",
                f"{len(multilayer_edges)} per-unit sharing edges exported, {cross} tagged crossing a known chromosomal cluster. This decomposes the aggregate sharing network by plasmid unit; it does not add new evidence.",
                ["network.multilayer_edges.tsv", "network.multilayer.graphml"], "Load the GraphML export to inspect which plasmid units drive cross-cluster sharing.", "sharing")

    samples = []
    for m in metadata:
        iso = m["isolate_id"]
        ps = [p for p in candidates if p["isolate_id"] == iso]
        ap = [p for p in ps if p.get("plasmid_unit")]
        calls = [f for f in eligible if f["isolate_id"] == iso]
        done = sum((p["plasmid_id"], "amr") in completed for p in ap)
        links = [e for e in payload["edges"] if iso in (e["source"], e["target"])]
        low = sum(p.get("quality_status") in {"low_confidence", "uncertain"} for p in ap)
        samples.append(dict(isolate_id=iso, organism=m.get("organism", ""), location=m.get("location", ""),
            date=m.get("date", ""), chromosomal_cluster=m.get("chromosomal_cluster", ""),
            submitted_candidates=len(ps), accepted_candidates=None if failed else len(ap),
            rejected_candidates=sum(p.get("quality_status") == "rejected" for p in ps),
            total_length_bp=None if failed else sum(int(p.get("length") or 0) for p in ap), units=None if failed else len({p["plasmid_unit"] for p in ap}),
            sharing_partners=None if failed else len(links), eligible_arg_labels=";".join(sorted({f["amr_gene"] for f in calls})),
            amr_completed=done, annotation_denominator=len(ap), quality_review_candidates=low,
            state="warn" if low or len(ps)>len(ap) else "info" if ap else "not_evaluated",
            interpretation=("Run incomplete: completed candidate assignments and sharing results are unavailable." if failed else
                           f"{len(ap)} accepted candidates; {len(links)} sharing partners; AMR stage complete for {done}/{len(ap)}. "
                            + ("Chromosome cluster available." if m.get("chromosomal_cluster") else "Chromosomal comparison not evaluated.")
                            + (" No indexed candidate does not demonstrate plasmid absence." if not ap else ""))))
    size = Counter()
    for p in accepted:
        length = int(p.get("length") or 0)
        size["<10 kb" if length<10000 else "10–<50 kb" if length<50000 else "50–<100 kb" if length<100000 else "≥100 kb"] += 1
    unit_members = Counter(p["plasmid_unit"] for p in accepted)
    functions = defaultdict(set)
    for f in features:
        functions[f.get("functional_category") or "unclassified"].add(f["plasmid_id"])
    drug = defaultdict(set)
    for f in eligible:
        drug[f.get("drug_class") or "Unspecified"].add(f["plasmid_id"])
    charts = dict(
        quality=[dict(label=k, value=v) for k,v in sorted(qc.items())],
        lengths=[dict(label=k, value=size[k]) for k in ("<10 kb","10–<50 kb","50–<100 kb","≥100 kb")],
        units=[dict(label=k, value=v) for k,v in sorted(unit_members.items(), key=lambda x:(-x[1],x[0]))],
        functions=[dict(label=k, value=len(v)) for k,v in sorted(functions.items())],
        drugs=[dict(label=k, value=len(v)) for k,v in sorted(drug.items())],
        metadata=[dict(label=k, value=len(metadata)-v, denominator=len(metadata)) for k,v in missing.items()],
        annotation=[dict(label=k, value=v, denominator=n) for k,v in completion.items()],
        timeline=[dict(label=k, value=v) for k,v in sorted(Counter(m.get("date","")[:7] for m in metadata if m.get("date")).items())],
        denominators=dict(submitted_candidates=total, accepted_candidates=n, isolates=len(metadata),
                          functions="Unique candidate carriers per category; categories overlap.",
                          drugs="Unique candidates carrying eligible ARGs per drug class; classes overlap.",
                          timeline="Dated isolates per collection month; not inferred transmission dates."))
    modules = [dict(module=f["title"], state=f["state"], interpretation=f["interpretation"],
                    evidence=";".join(f["evidence"]), next_step=f["next_step"], section=f["section"]) for f in findings]
    return dict(rule_version=RULE_VERSION, findings=findings, samples=samples, modules=modules, charts=charts)
