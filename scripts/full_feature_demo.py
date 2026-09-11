#!/usr/bin/env python3
"""Full-feature synthetic demo: exercises every dashboard capability the minimal seeded
`plasmicord demo` skips in one run -- external typing, containment, size-corrected
clustering, a multilayer network export, and quality-evidence copy_number -- on
illustrative seeded random sequences, not biological detections.

Cohort design (8 isolates, 4 chromosomal clusters):
  iso1 (CC1): A, D          iso5 (CC3): no plasmids, dated
  iso2 (CC1): A, B          iso6 (CC4): C_full, undated
  iso3 (CC2): A' (~0.5% mutated A -> same unit as A: CC1<->CC2 cross-cluster link)
  iso4 (CC3): B' (~0.5% mutated B -> same unit as B: CC1<->CC3 cross-cluster link)
  iso7 (CC4): C_short (75% length prefix of C_full -> containment candidate, missing location)
  iso8 (CC2): D_longer (60% longer than D -> merges with D only under size correction:
              CC1<->CC2 link visible only with --size-correction-per-percent, missing organism)

Parameters (k=13, threshold=0.01, containment-max-distance=0.05, size-correction-per-percent=0.01,
size-correction-cap-pct=40) were chosen empirically against this project's own k-mer distance
formula (python/kmer_distance.py) so every optional dataset in the report actually populates;
see the plan/verification notes in the accompanying commit for the tuning values used.
"""
import argparse
from pathlib import Path
import random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from python.cli import run
from python.contracts import write_tsv, FEATURE_FIELDS
from python.population_summary import summarize


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, rate, rng):
    return "".join(rng.choice([c for c in "ACGT" if c != base]) if rng.random() < rate else base for base in seq)


def digest(seq):
    import hashlib
    return hashlib.sha256(seq.encode()).hexdigest()


def build_sequences(seed=7):
    rng = random.Random(seed)
    A = rand_seq(1200, rng)
    Ap = mutate(A, 0.005, rng)
    B = rand_seq(900, rng)
    Bp = mutate(B, 0.005, rng)
    C_full = rand_seq(800, rng)
    C_short = C_full[: int(len(C_full) * 0.75)]
    D = rand_seq(500, rng)
    D_longer = D + rand_seq(int(len(D) * 0.6), rng)
    return dict(A=A, Ap=Ap, B=B, Bp=Bp, C_full=C_full, C_short=C_short, D=D, D_longer=D_longer)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=Path("results_full_demo"))
    ap.add_argument("--engine", choices=("kmer", "mash"), default="kmer")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    base = args.out.resolve()
    if base.exists() and any(base.iterdir()):
        ap.error(f"Output directory is not empty: {base}. Use a fresh directory.")
    inputs = base / "inputs"
    inputs.mkdir(parents=True)
    seq = build_sequences(args.seed)

    plasmids = [
        ("iso1", "A", seq["A"]), ("iso1", "D", seq["D"]),
        ("iso2", "A", seq["A"]), ("iso2", "B", seq["B"]),
        ("iso3", "Ap", seq["Ap"]),
        ("iso4", "Bp", seq["Bp"]),
        ("iso6", "Cfull", seq["C_full"]),
        ("iso7", "Cshort", seq["C_short"]),
        ("iso8", "Dlonger", seq["D_longer"]),
    ]
    manifest_rows = []
    for iso, name, sequence in plasmids:
        pid = f"{iso}__{name}"
        (inputs / f"{pid}.fasta").write_text(f">{pid}\n{sequence}\n", encoding="utf-8")
        manifest_rows.append(dict(isolate_id=iso, plasmid_id=pid, fasta_path=f"{pid}.fasta",
                                  source_tool="synthetic_fixture", source_tool_version="1"))
    write_tsv(inputs / "manifest.tsv", ["isolate_id", "plasmid_id", "fasta_path", "source_tool", "source_tool_version"],
              manifest_rows)

    metadata_rows = [
        dict(isolate_id="iso1", chromosomal_cluster="CC1", location="Site A", date="2026-01-01", organism="Synthetic organism"),
        dict(isolate_id="iso2", chromosomal_cluster="CC1", location="Site A", date="2026-01-15", organism="Synthetic organism"),
        dict(isolate_id="iso3", chromosomal_cluster="CC2", location="Site B", date="2026-02-01", organism="Synthetic organism"),
        dict(isolate_id="iso4", chromosomal_cluster="CC3", location="Site C", date="2026-02-20", organism="Synthetic organism"),
        dict(isolate_id="iso5", chromosomal_cluster="CC3", location="Site C", date="2026-03-01", organism="Synthetic organism"),
        dict(isolate_id="iso6", chromosomal_cluster="CC4", location="Site D", date="", organism="Synthetic organism"),
        dict(isolate_id="iso7", chromosomal_cluster="CC4", location="", date="2026-03-10", organism="Synthetic organism"),
        dict(isolate_id="iso8", chromosomal_cluster="CC2", location="Site B", date="2026-03-20", organism=""),
    ]
    write_tsv(inputs / "metadata.tsv", ["isolate_id", "chromosomal_cluster", "location", "date", "organism"], metadata_rows)

    features_rows = [
        dict(plasmid_id="iso1__A", feature_id="iso1__A_bla", gene_symbol="blaTEM-1", product_name="beta-lactamase TEM-1 (illustrative)",
             start=40, end=780, strand="+", functional_category="amr", amr_gene="blaTEM-1", drug_class="beta-lactam",
             hit_class="curated", annotation_confidence="high", annotation_engine="synthetic_fixture",
             database_name="synthetic_not_a_biological_database", database_version="1"),
        dict(plasmid_id="iso2__B", feature_id="iso2__B_mer", gene_symbol="merA", product_name="mercuric reductase (illustrative)",
             start=30, end=420, strand="-", functional_category="metal_resistance"),
        dict(plasmid_id="iso6__Cfull", feature_id="iso6__Cfull_is", gene_symbol="IS26", product_name="insertion sequence IS26 (illustrative)",
             start=50, end=350, strand="+", functional_category="insertion_sequence"),
        dict(plasmid_id="iso1__D", feature_id="iso1__D_mge", gene_symbol="tnpA", product_name="transposase (illustrative)",
             start=20, end=250, strand="+", functional_category="mobile_genetic_element"),
        dict(plasmid_id="iso2__A", feature_id="iso2__A_rep", gene_symbol="repA", product_name="replication initiator (illustrative)",
             start=800, end=1100, strand="+", functional_category="replication"),
    ]
    write_tsv(inputs / "features.tsv", FEATURE_FIELDS, features_rows)

    quality_evidence_rows = [
        dict(plasmid_id="iso1__A", sequence_sha256=digest(seq["A"]), evidence_source="illustrative_read_mapping",
             classification="plasmid", read_breadth="0.98", mean_depth="25", chromosome_fraction="0.01", copy_number="3.5"),
    ]
    write_tsv(inputs / "quality_evidence.tsv",
              ["plasmid_id", "sequence_sha256", "evidence_source", "classification", "read_breadth",
               "mean_depth", "chromosome_fraction", "copy_number"], quality_evidence_rows)

    external_typing_rows = [
        dict(plasmid_id="iso1__A", sequence_sha256=digest(seq["A"]), evidence_source="illustrative_mob_suite_run",
             external_tool="mob_suite", mob_primary_cluster_id="AA123", ptu_assignment="PTU-FE",
             ptu_confidence="high", predicted_host_range_overall_name="Enterobacteriaceae",
             associated_pmids="12345678", predicted_transmissibility_score="0.82",
             predicted_transmissibility_call="conjugative", plasmidfinder_inc_types="IncF",
             plasmidfinder_identity="98.5", predicted_classification_call="plasmid",
             predicted_classification_score="0.95", plsdb_nearest_accession="NZ_CP000000.1",
             plsdb_nearest_distance="0.02", plsdb_nearest_host="Escherichia coli"),
    ]
    from python.external_typing import EXTERNAL_TYPING_FIELDS
    write_tsv(inputs / "external_typing.tsv", EXTERNAL_TYPING_FIELDS, external_typing_rows)

    out = base / "results"
    run(argparse.Namespace(manifest=inputs / "manifest.tsv", metadata=inputs / "metadata.tsv",
                           features=inputs / "features.tsv", annotation_config=None, annotation_profile="essential",
                           annotation_cache=None, threads=1,
                           quality_evidence=inputs / "quality_evidence.tsv",
                           external_typing=inputs / "external_typing.tsv",
                           out=out, mode="precomputed", engine=args.engine, threshold=0.01, linkage="single",
                           k=13, sketch_size=1000, min_length=50,
                           containment_min_ratio=0.5, containment_max_ratio=0.95, containment_max_distance=0.05,
                           multilayer_network=True, size_correction_per_percent=0.01, size_correction_cap_pct=40.0,
                           synthetic=True))
    summarize(argparse.Namespace(results=out, out_prefix=None))
    print(f"Full-feature demo complete: {out / 'REPORT.html'}")
    print("Illustrative synthetic sequences and gene annotations; not biological detections.")


if __name__ == "__main__":
    main()
