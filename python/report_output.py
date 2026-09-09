"""Report rendering, portable bundles and regeneration without rerunning analysis."""
from datetime import datetime, timezone
import html
import json
import re
from pathlib import Path
import shutil
import zipfile
from .contracts import read_tsv, write_tsv
from .report_evidence import RULE_VERSION, PRESENTATION, catalog, describe, result_files, sha256, summary
from . import __version__

GENERATED = PRESENTATION | {"interpretations.json", "sample_summary.tsv", "module_summary.tsv", "report_charts.json"}


def attach_calibration(out, source, provenance):
    if not source:
        source = out / "calibration"
    source = Path(source).resolve()
    if not source.is_dir():
        return None
    record_path = source/"calibration.json"
    if not record_path.is_file():
        raise ValueError("Calibration attachment requires calibration.json")
    record = json.loads(record_path.read_text(encoding="utf-8"))
    matrix = out/"plasmid_matrix.tsv"
    if not matrix.is_file() or record.get("matrix_sha256") != sha256(matrix):
        raise ValueError("Calibration matrix does not match this run; refusing to mix cohorts")
    if record.get("distance_engine") != provenance.get("distance_engine") or record.get("linkage") != provenance.get("linkage"):
        raise ValueError("Calibration engine/linkage does not match this run")
    for key, value in record.get("source_parameters", {}).items():
        if key in provenance and provenance[key] != value:
            raise ValueError(f"Calibration parameter mismatch: {key}")
    required = ("calibration.json", "calibration_metrics.tsv", "calibration_labels.tsv")
    for name in required:
        if not (source/name).is_file():
            raise ValueError(f"Calibration attachment missing {name}")
    destination = out/"calibration"
    if source != destination.resolve():
        if destination.exists() and any(destination.iterdir()):
            raise ValueError("A calibration attachment already exists; preserve it or use a separate result copy")
        destination.mkdir(exist_ok=True)
        for name in (*required, "CALIBRATION.html", "CALIBRATION.md"):
            if (source/name).is_file():
                shutil.copy2(source/name, destination/name)
    return dict(record=record, metrics=read_tsv(destination/"calibration_metrics.tsv"))


def render_dashboard(out, payload, calibration_source=None, bundle=True):
    out = Path(out).resolve()
    provenance = payload["provenance"]
    if not bundle and (out/"REPORT_BUNDLE.zip").is_file():
        (out/"REPORT_BUNDLE.zip").unlink()
    payload["calibration"] = attach_calibration(out, calibration_source, provenance)
    generated_at = datetime.now(timezone.utc).isoformat()
    payload["report"] = dict(schema_version="2.0", renderer_version=__version__, generated_at=generated_at, rule_version=RULE_VERSION,
        scope="Cohort statistics describe the full available run; selection-specific interpretations are explicitly labelled.",
        snapshot_note="Embedded file previews and metadata are bounded snapshots. output_manifest.json records final file checksums; raw files remain authoritative.",
        bundle=bundle, offline=True)
    payload["dashboard"] = summary(payload)
    dashboard = payload["dashboard"]
    (out/"interpretations.json").write_text(json.dumps(dict(rule_version=RULE_VERSION, findings=dashboard["findings"]), indent=2), encoding="utf-8")
    (out/"report_charts.json").write_text(json.dumps(dashboard["charts"], indent=2), encoding="utf-8")
    write_tsv(out/"sample_summary.tsv", list(dashboard["samples"][0]) if dashboard["samples"] else ["isolate_id","state","interpretation"], dashboard["samples"])
    write_tsv(out/"module_summary.tsv", ["module","state","interpretation","evidence","next_step","section"], dashboard["modules"])
    # The in-progress run record must not appear as a stale "running" preview.
    provenance["report_generated_at"] = generated_at
    provenance["report_rule_version"] = RULE_VERSION
    provenance["report_bundle"] = bundle
    (out/"run_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    files = catalog(out)
    # These products contain the inventory, so their final sizes/hashes cannot be
    # embedded recursively. They are explicit deferred entries, never fake hashes.
    for name in sorted(PRESENTATION):
        if name == "REPORT_BUNDLE.zip" and not bundle:
            continue
        module, description = describe(name)
        files.append(dict(path=name, href="./"+name, name=name, module=module, description=description,
                          size_bytes=None, sha256=None, modified_utc=generated_at, media_type=Path(name).suffix[1:],
                          preview=None, deferred=True))
    payload["artifacts"] = sorted(files, key=lambda f:f["path"].lower())
    # Embed a bounded heatmap; the full matrix is always downloadable.
    matrix = out/"plasmid_matrix.tsv"
    payload["distance_view"] = dict(ids=[], matrix=[], limited=False, limit=120)
    if matrix.is_file() and provenance.get("status") == "complete":
        from .cluster_plasmids import read_matrix
        ids, _, values = read_matrix(matrix)
        payload["distance_view"] = dict(ids=ids[:120], matrix=[row[:120] for row in values[:120]],
                                       limited=len(ids)>120, total=len(ids), limit=120)
    data = json.dumps(payload).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    links = "".join(f'<li><a download href="{html.escape(f["href"], quote=True)}">{html.escape(f["path"])}</a></li>' for f in payload["artifacts"])
    folder = Path(__file__).parent
    template = (folder/"report_template.html").read_text(encoding="utf-8")
    # One-pass substitution keeps input strings from becoming template directives.
    javascript = (folder/"report.js").read_text(encoding="utf-8")+"\n"+(folder/"report_dashboard.js").read_text(encoding="utf-8")
    replacements = {"REPORT_CSS":(folder/"report.css").read_text(encoding="utf-8"),
                    "REPORT_JS":javascript, "DOWNLOADS":links, "REPORT_DATA":data}
    page = re.sub(r"__(REPORT_CSS|REPORT_JS|DOWNLOADS|REPORT_DATA)__",
                  lambda match: replacements[match.group(1)], template)
    (out/"report_data.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (out/"REPORT.html").write_text(page, encoding="utf-8")
    lines = ["# "+provenance["project_title"], "", provenance["project_tagline"], "",
             f"Report generated: {generated_at}. Rule set: {RULE_VERSION}. Run state: {provenance.get('status')}.",
             "", "## Cohort counts", "", *[f"- {k}: {v if v is not None else 'not evaluated'}" for k,v in payload["counts"].items()],
             "", "## Automated interpretation", ""]
    for f in dashboard["findings"]:
        lines.extend([f"### {f['title']} — {f['state']}", "", f["interpretation"], "",
                      "Review: "+f["next_step"], "", "Evidence: "+", ".join(f["evidence"]), ""])
    lines.extend(["## Individual isolates", ""])
    for sample in dashboard["samples"]:
        lines.extend([f"### {sample['isolate_id']}", "", sample["interpretation"], ""])
    lines.extend(["## Methods and provenance", "", "Run parameters and tool versions:", "", json.dumps(provenance, indent=2),
                  "", "## Interpretation limits", "", *["- "+s for s in payload["safeguards"]],
                  "## Output inventory", "", *[f"- [{f['path']}]({f['href']}): {f['description']}" for f in payload["artifacts"]]])
    (out/"REPORT.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    stable = [p for p in result_files(out) if p.relative_to(out).as_posix() not in {"run_provenance.json","output_manifest.json","REPORT_BUNDLE.zip"}]
    provenance["output_sha256"] = {p.relative_to(out).as_posix():sha256(p) for p in stable}
    (out/"run_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    inventory = catalog(out, include_generated=True, previews=False)
    manifest = dict(schema_version="1.0", generated_at=generated_at, files=inventory,
                    excluded_from_hashes=["output_manifest.json","REPORT_BUNDLE.zip"],
                    exclusion_reason="Manifest and bundle contain this inventory; self-referential hashes are not defined.",
                    file_scope="Regular result-directory files only; symlinks and paths outside the result directory are excluded.")
    (out/"output_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if bundle:
        archive = out/"REPORT_BUNDLE.zip"
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, strict_timestamps=False) as z:
            for p in result_files(out):
                if p != archive:
                    z.write(p, p.relative_to(out).as_posix())
    return payload


def refresh_report(args):
    out = args.results.resolve()
    run = json.loads((out/"run_provenance.json").read_text(encoding="utf-8"))
    payload = json.loads((out/"report_data.json").read_text(encoding="utf-8"))
    # Regeneration changes presentation, not the meaning of an altered data set.
    for rel, digest in run.get("output_sha256", {}).items():
        if rel in GENERATED:
            continue
        path = out/rel
        if not path.resolve().is_relative_to(out) or not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Result changed or missing since the recorded run: {rel}")
    expected = run.get("output_sha256", {}).get("report_data.json")
    if expected and sha256(out/"report_data.json") != expected:
        raise ValueError("Report data checksum differs from the recorded run")
    payload["provenance"] = run
    if getattr(args, "calibration", None) and not args.calibration.is_dir():
        raise ValueError("Calibration attachment directory does not exist")
    render_dashboard(out, payload, getattr(args, "calibration", None), not getattr(args, "no_bundle", False))
    print(f"Report regenerated without rerunning analysis: {out/'REPORT.html'}")


def failure_report(out, provenance, index, meta, sequences):
    from .report import SAFEGUARDS
    payload = dict(schema_version="1.0", provenance=provenance, metadata=meta,
        plasmids=[dict(p, plasmid_unit="") for p in index], features=[], unit_functions=[], edges=[],
        edge_evidence=[], crosslinks=[], safeguards=SAFEGUARDS,
        biological_quality=read_tsv(out/"biological_quality.tsv") if (out/"biological_quality.tsv").is_file() else [],
        annotation_status=read_tsv(out/"annotation_status.tsv") if (out/"annotation_status.tsv").is_file() else [],
        plasmid_clusters=read_tsv(out/"plasmid_clusters.tsv") if (out/"plasmid_clusters.tsv").is_file() else [],
        contigs={pid:[dict(id=n,length=len(s)) for n,s in records] for pid,records in sequences.items()},
        sensitivity=[], counts=dict(isolates=len(meta),plasmids=None,units=None,sharing_pairs=None,
                                   cross_cluster_pairs=None,cross_cluster_unit_links=None,typed_isolates=sum(bool(m.get("chromosomal_cluster")) for m in meta),
                                   observed_arg_plasmids=None))
    # Partial distances may be malformed. Do not attempt a completed-run heatmap.
    return render_dashboard(out, payload, bundle=False)


def add_parser(subs):
    p = subs.add_parser("report", help="Regenerate a detailed offline dashboard from a recorded result directory")
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--calibration", type=Path, help="Attach a matching calibration directory; no thresholds are changed")
    p.add_argument("--no-bundle", action="store_true", help="Do not create a result ZIP, useful for very large runs")
    p.set_defaults(func=refresh_report)
