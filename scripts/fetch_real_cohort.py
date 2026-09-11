#!/usr/bin/env python3
"""Download a small, real, accession-backed plasmid cohort from NCBI (stdlib urllib only,
no third-party HTTP client) and lay it out exactly as scripts/prepare_reference_validation.py
expects, so that existing script can build a manifest unchanged.

This is the only network-touching code in the repository. It is deliberately kept out of
`plasmicord run`/`demo` and out of the default test suite: PlasmiCord's core analysis requires
no network access. Use it only when you explicitly want a live pipeline-completeness check
against real public sequences (see test/test_real_cohort.py and `make real-cohort-test`).

Default cohort: the three smallest BioProject PRJNA636382 (Snyder et al. 2020, recurrent-UTI
E. coli) sister isolates already vetted in docs/validation/reference_metadata.tsv -- about
276 KB of plasmid sequence across 6 plasmids total. Override with --assemblies.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from python.contracts import checksum, write_tsv

DEFAULT_ASSEMBLIES = {
    "upec_ecpf5": "GCF_013372425.1",
    "upec_ecpf7": "GCF_013372405.1",
    "upec_ecpf14": "GCF_013372385.1",
}
ORGANISM = "Escherichia coli"
BIOPROJECT = "PRJNA636382"
SOURCE_STUDY = "Snyder_2020_recurrent_UTI"

DATASETS_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/sequence_reports"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REQUEST_DELAY_SECONDS = 0.4  # a fixed courtesy delay between requests; well under NCBI's unauthenticated rate limit


def fetch_json(url, contact_email=None):
    query = dict(tool="plasmicord")
    if contact_email:
        query["email"] = contact_email
    full_url = url + ("&" if "?" in url else "?") + urllib.parse.urlencode(query)
    request = urllib.request.Request(full_url, headers={"User-Agent": "plasmicord-fetch-real-cohort (stdlib urllib)"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError) as error:
        raise ValueError(f"Could not reach NCBI at {url}: {error}") from error


def fetch_text(url, params, contact_email=None):
    query = dict(params, tool="plasmicord")
    if contact_email:
        query["email"] = contact_email
    full_url = url + "?" + urllib.parse.urlencode(query)
    request = urllib.request.Request(full_url, headers={"User-Agent": "plasmicord-fetch-real-cohort (stdlib urllib)"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return response.read().decode("utf-8")
    except (urllib.error.URLError, urllib.error.HTTPError) as error:
        raise ValueError(f"Could not reach NCBI at {url}: {error}") from error


def plasmid_reports(assembly_accession, contact_email=None):
    payload = fetch_json(DATASETS_URL.format(accession=assembly_accession), contact_email)
    reports = payload.get("reports", [])
    if not reports:
        raise ValueError(f"NCBI Datasets returned no sequence reports for {assembly_accession}")
    return [r for r in reports if r.get("assigned_molecule_location_type") == "Plasmid"]


def fetch_cohort(assemblies, out, contact_email=None):
    out = Path(out).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Output directory is not empty: {out}. Use a fresh directory.")
    references = out / "references"
    references.mkdir(parents=True)
    cohort_rows, download_records = [], []
    for sample_id, assembly_accession in assemblies.items():
        reports = plasmid_reports(assembly_accession, contact_email)
        time.sleep(REQUEST_DELAY_SECONDS)
        if not reports:
            raise ValueError(f"{assembly_accession}: NCBI reports no plasmid molecules for this assembly")
        sample_dir = references / sample_id
        sample_dir.mkdir(parents=True)
        fasta_path = sample_dir / "reference.fna"
        sequence_report_path = sample_dir / "sequence_report.jsonl"
        with fasta_path.open("w", encoding="utf-8") as fasta_handle, \
             sequence_report_path.open("w", encoding="utf-8") as report_handle:
            for report in reports:
                accession = report.get("refseq_accession") or report.get("genbank_accession")
                fasta = fetch_text(EFETCH_URL, dict(db="nuccore", id=accession, rettype="fasta", retmode="text"), contact_email)
                time.sleep(REQUEST_DELAY_SECONDS)
                if not fasta.startswith(">"):
                    raise ValueError(f"{accession}: efetch did not return FASTA text")
                fasta_handle.write(fasta if fasta.endswith("\n") else fasta + "\n")
                # Translate the REST API's snake_case fields to the camelCase fields
                # scripts/prepare_reference_validation.py already expects (NCBI Datasets
                # CLI's own JSONL export convention) so that script can be reused unchanged.
                report_handle.write(json.dumps(dict(
                    refseqAccession=report.get("refseq_accession"),
                    genbankAccession=report.get("genbank_accession"),
                    assignedMoleculeLocationType=report.get("assigned_molecule_location_type"),
                    assemblyAccession=report.get("assembly_accession"),
                    length=report.get("length"),
                )) + "\n")
                download_records.append(dict(sample_id=sample_id, accession=accession,
                                             source_url=EFETCH_URL, length=report.get("length")))
        cohort_rows.append(dict(sample_id=sample_id, assembly_accession=assembly_accession,
                                organism=ORGANISM, bioproject=BIOPROJECT, source_study=SOURCE_STUDY,
                                truth_technology=""))
        download_records.append(dict(sample_id=sample_id, accession=assembly_accession,
                                     source_url=DATASETS_URL.format(accession=assembly_accession), length=None))
    write_tsv(out / "cohort.tsv", ["sample_id", "assembly_accession", "organism", "bioproject",
                                  "source_study", "truth_technology"], cohort_rows)
    manifest = dict(schema_version="1.0", source="NCBI Datasets v2 REST API + E-utilities efetch",
                    assemblies=assemblies, downloads=download_records,
                    reference_sha256={str(p.relative_to(out)): checksum(p) for p in references.rglob("*") if p.is_file()},
                    limitation="Real public reference sequences; no epidemiological cluster, date "
                               "or location is invented from sequence similarity.")
    (out / "download_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Downloaded {len(cohort_rows)} isolates to {out}")
    return out


def parse_assemblies(text):
    assemblies = {}
    for entry in text.split(","):
        sample_id, _, accession = entry.partition(":")
        if not sample_id or not accession:
            raise ValueError(f"Invalid --assemblies entry (expected sample_id:accession): {entry}")
        assemblies[sample_id] = accession
    return assemblies


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--assemblies", help="Comma-separated sample_id:accession pairs; default is a small "
                                        "already-vetted 3-isolate BioProject PRJNA636382 subset")
    ap.add_argument("--contact-email", help="Optional contact email sent to NCBI per E-utilities usage guidelines")
    args = ap.parse_args()
    assemblies = parse_assemblies(args.assemblies) if args.assemblies else DEFAULT_ASSEMBLIES
    try:
        fetch_cohort(assemblies, args.out, args.contact_email)
    except ValueError as error:
        ap.exit(1, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
