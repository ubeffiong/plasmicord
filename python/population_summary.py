"""Optional post-processing: cohort-level aggregation over an already-completed run.

Pure aggregation over outputs `plasmicord run` already wrote; no new upstream data and
no new inference. Adapted to PlasmiCord's own metadata contract: the free-text `location`
column is used in place of a dedicated country field (phac-nml/plasmid_analysis's schema).
"""
import json
from datetime import date
from pathlib import Path
from .contracts import checksum, read_tsv, write_tsv

PU_LEVEL_FIELDS = ("plasmid_unit n_plasmids n_isolates isolates replicon_types mobility_classes "
                   "n_organisms organisms n_locations locations first_date last_date date_range_days "
                   "n_resistance_genes resistance_genes n_drug_classes drug_classes interpretation").split()

DIMENSION_FIELDS = ("dimension value n_isolates n_plasmids n_plasmid_units "
                    "n_resistance_genes resistance_genes first_date last_date interpretation").split()

SOURCE_FILES = ("plasmid_index.tsv", "metadata.tsv", "plasmid_clusters.tsv",
                "biological_quality.tsv", "functional_features.tsv")

PU_INTERPRETATION = "aggregation of already-computed run outputs; not a validated core/accessory or transmission claim"
DIMENSION_INTERPRETATION = "aggregation of already-computed run outputs across a metadata dimension; not a validated epidemiological claim"
DIMENSIONS = ("organism", "location", "chromosomal_cluster")


def _date_range(dates):
    dates = sorted(d for d in dates if d)
    if not dates:
        return "", "", ""
    first, last = dates[0], dates[-1]
    try:
        days = (date.fromisoformat(last) - date.fromisoformat(first)).days
    except ValueError:
        days = ""
    return first, last, days


def _headline_genes(features, plasmid_ids, field):
    return sorted({f[field] for f in features
                  if f["plasmid_id"] in plasmid_ids and f.get("headline_eligible") == "true" and f.get(field)})


def summarize(args):
    out = args.results.resolve()
    provenance_path = out / "run_provenance.json"
    if not provenance_path.is_file():
        raise ValueError(f"{out}: not a PlasmiCord result directory (missing run_provenance.json)")
    run = json.loads(provenance_path.read_text(encoding="utf-8"))
    if run.get("status") != "complete":
        raise ValueError("population-summary requires a completed run")
    output_sha256 = run.get("output_sha256", {})
    for rel in SOURCE_FILES:
        path = out / rel
        expected = output_sha256.get(rel)
        if not path.is_file() or (expected and checksum(path) != expected):
            raise ValueError(f"Result changed or missing since the recorded run: {rel}")

    index = read_tsv(out / "plasmid_index.tsv")
    meta_by_iso = {r["isolate_id"]: r for r in read_tsv(out / "metadata.tsv")}
    clusters = {r["plasmid_id"]: r["plasmid_unit"] for r in read_tsv(out / "plasmid_clusters.tsv")}
    quality_by_pid = {r["plasmid_id"]: r for r in read_tsv(out / "biological_quality.tsv")}
    features = read_tsv(out / "functional_features.tsv")

    pu_plasmids = {}
    for row in index:
        pu = clusters.get(row["plasmid_id"])
        if pu:
            pu_plasmids.setdefault(pu, []).append(row)

    prefix = Path(getattr(args, "out_prefix", None) or "population_summary")
    if not prefix.is_absolute():
        prefix = out / prefix

    pu_rows = []
    for pu, plasmids in sorted(pu_plasmids.items()):
        pids = {p["plasmid_id"] for p in plasmids}
        isolates = sorted({p["isolate_id"] for p in plasmids})
        organisms = sorted({meta_by_iso[i]["organism"] for i in isolates if meta_by_iso.get(i, {}).get("organism")})
        locations = sorted({meta_by_iso[i]["location"] for i in isolates if meta_by_iso.get(i, {}).get("location")})
        replicons = sorted({quality_by_pid[pid]["replicon_type"] for pid in pids
                            if quality_by_pid.get(pid, {}).get("replicon_type")})
        mobility = sorted({quality_by_pid[pid]["mobility_class"] for pid in pids
                           if quality_by_pid.get(pid, {}).get("mobility_class") not in (None, "", "not_evaluated")})
        first_date, last_date, days = _date_range(meta_by_iso.get(i, {}).get("date", "") for i in isolates)
        resistance = _headline_genes(features, pids, "amr_gene")
        drug_classes = _headline_genes(features, pids, "drug_class")
        pu_rows.append(dict(plasmid_unit=pu, n_plasmids=len(pids), n_isolates=len(isolates),
                            isolates=";".join(isolates), replicon_types=";".join(replicons),
                            mobility_classes=";".join(mobility), n_organisms=len(organisms),
                            organisms=";".join(organisms), n_locations=len(locations), locations=";".join(locations),
                            first_date=first_date, last_date=last_date, date_range_days=days,
                            n_resistance_genes=len(resistance), resistance_genes=";".join(resistance),
                            n_drug_classes=len(drug_classes), drug_classes=";".join(drug_classes),
                            interpretation=PU_INTERPRETATION))
    write_tsv(Path(f"{prefix}.pu_level.tsv"), PU_LEVEL_FIELDS, pu_rows)

    dimension_rows = []
    for dimension in DIMENSIONS:
        values = sorted({r.get(dimension) for r in meta_by_iso.values() if r.get(dimension)})
        for value in values:
            isolates = sorted(i for i, r in meta_by_iso.items() if r.get(dimension) == value)
            plasmids = [p for p in index if p["isolate_id"] in isolates]
            pids = {p["plasmid_id"] for p in plasmids}
            units = sorted({clusters[p["plasmid_id"]] for p in plasmids if clusters.get(p["plasmid_id"])})
            first_date, last_date, _ = _date_range(meta_by_iso[i].get("date", "") for i in isolates)
            resistance = _headline_genes(features, pids, "amr_gene")
            dimension_rows.append(dict(dimension=dimension, value=value, n_isolates=len(isolates),
                                       n_plasmids=len(pids), n_plasmid_units=len(units),
                                       n_resistance_genes=len(resistance), resistance_genes=";".join(resistance),
                                       first_date=first_date, last_date=last_date, interpretation=DIMENSION_INTERPRETATION))
    write_tsv(Path(f"{prefix}.metadata_dimension.tsv"), DIMENSION_FIELDS, dimension_rows)
    print(f"Population summary: {prefix}.pu_level.tsv, {prefix}.metadata_dimension.tsv")


def add_parser(subs):
    p = subs.add_parser("population-summary", help="Aggregate a completed run's outputs into cohort-level summary tables")
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--out-prefix", help="Output path prefix (default: <results>/population_summary)")
    p.set_defaults(func=summarize)
