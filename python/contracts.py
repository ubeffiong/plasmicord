"""Input validation and normalized provenance contracts (stdlib only)."""
import csv
import hashlib
import re
from datetime import date
from pathlib import Path

MODES = ("precomputed", "longread", "hybrid", "plasbench")
CIRCULARITY = {"confirmed", "assembly_supported", "tool_reported", "linear", "unresolved"}
QUALITY = {"high_confidence", "moderate_confidence", "low_confidence", "uncertain", "rejected"}
INDEX_FIELDS = "isolate_id plasmid_id fasta_path source_type source_tool source_tool_version source_sequence_id sequencing_technology assembly_method polishing_method circularity_status declared_quality_status quality_status quality_warnings sequence_sha256 length n_contigs duplicate_of annotation_path reads_path assembly_graph_path".split()
FEATURE_FIELDS = "isolate_id plasmid_id plasmid_unit feature_id gene_symbol product_name feature_type start end strand functional_category functional_subcategory amr_gene drug_class resistance_mechanism replicon_type mobility_function ko_id kegg_module cog_category go_terms ec_number pfam_ids identity coverage hit_class annotation_confidence annotation_engine annotation_engine_version database_name database_version sequence_sha256 source_sequence_id headline_eligible detection_method reference_accession dbxref annotation_scope".split()


def read_tsv(path, required=()):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader((line for line in handle if line.strip() and not line.startswith("#")), delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or any(c not in fields for c in required):
            raise ValueError(f"{path}: missing required or duplicate columns; required: {', '.join(required)}")
        rows = []
        for row in reader:
            if None in row:
                raise ValueError(f"{path}: too many columns in a row")
            rows.append({k: (v or "").strip() for k, v in row.items()})
        return rows


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fasta_records(path):
    records, seen, name, parts = [], set(), None, []
    def save():
        if name is not None:
            seq = "".join(parts).upper()
            if not seq or set(seq) - set("ACGTRYSWKMBDHVN"):
                raise ValueError(f"{path}: empty sequence or invalid nucleotide alphabet in {name}")
            records.append((name, seq))
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                tokens = line[1:].split()
                if not tokens or tokens[0] in seen:
                    raise ValueError(f"{path}: empty or duplicate FASTA identifier")
                name, parts = tokens[0], []
                seen.add(name)
            elif name is None:
                raise ValueError(f"{path}: expected FASTA header; raw FASTQ is not a plasmid assembly")
            else:
                parts.append(line)
    save()
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    return records


def metadata(path):
    rows = read_tsv(path, ("isolate_id",))
    seen = set()
    for row in rows:
        iso = row["isolate_id"]
        if not iso or iso in seen:
            raise ValueError("Metadata isolate IDs must be unique and non-empty")
        seen.add(iso)
        if row.get("date"):
            date.fromisoformat(row["date"])
        if row.get("chromosomal_cluster", "").upper() in {"UNKNOWN", "NA", "N/A", "NONE"}:
            row["chromosomal_cluster"] = ""
    if not rows:
        raise ValueError("Metadata must contain at least one isolate")
    return rows


def validate_manifest(path, meta, mode, min_length=200):
    rows = read_tsv(path, ("isolate_id", "plasmid_id", "fasta_path"))
    seen, hashes, sequences = set(), {}, {}
    isolates = {r["isolate_id"] for r in meta}
    for row in rows:
        pid = row["plasmid_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", pid) or pid.casefold() in seen or pid.endswith(".") or pid.upper().split(".")[0] in {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)), *(f"LPT{i}" for i in range(1, 10))}:
            raise ValueError(f"Invalid, reserved or duplicate plasmid_id: {pid}")
        seen.add(pid.casefold())
        if row["isolate_id"] not in isolates:
            raise ValueError(f"{pid}: isolate is missing from metadata")
        fasta = Path(row["fasta_path"])
        if not fasta.is_absolute():
            fasta = Path(path).resolve().parent / fasta
        records = fasta_records(fasta)
        for optional_path in ('assembly_graph_path', 'annotation_path', 'reads_path'):
            if row.get(optional_path):
                p = Path(row[optional_path])
                row[optional_path] = str(p if p.is_absolute() else Path(path).resolve().parent / p)
        sequences[pid] = records
        # Sequence-only digest preserves record boundaries and ignores identifiers/wrapping.
        digest = hashlib.sha256("\n".join(seq for _, seq in records).encode()).hexdigest()
        if row.get("sequence_sha256") and row["sequence_sha256"] != digest:
            raise ValueError(f"{pid}: sequence checksum mismatch")
        length = sum(len(seq) for _, seq in records)
        circularity = row.get("circularity_status") or "unresolved"
        declared = row.get("quality_status") or "uncertain"
        if circularity not in CIRCULARITY or declared not in QUALITY:
            raise ValueError(f"{pid}: invalid circularity or quality status")
        warnings = ["plasmid_identity_not_independently_validated"]
        status = "uncertain"  # do not promote a self-declared quality label
        if len(records) > 1:
            warnings.append("fragmented_reconstruction")
            status = "low_confidence"
        if circularity == "confirmed":
            warnings.append("circularity_confirmation_is_submitter_supplied")
        tech = row.get("sequencing_technology", "").lower()
        if ("ont" in tech or "nanopore" in tech) and row.get("polishing_method", "").lower() in {"", "none", "unknown"}:
            warnings.append("ont_polishing_not_documented")
            status = "low_confidence"
        if any(set(seq) - set("ACGT") for _, seq in records):
            warnings.append("ambiguous_bases")
            status = "low_confidence"
        if length < min_length or declared == "rejected":
            status = "rejected"
            warnings.append("below_minimum_length" if length < min_length else "submitter_rejected")
        if declared == "low_confidence" and status != "rejected":
            status = "low_confidence"
        row.update(fasta_path=str(fasta.resolve()), source_type=row.get("source_type") or mode,
                   circularity_status=circularity, declared_quality_status=declared,
                   quality_status=status, quality_warnings=";".join(warnings),
                   source_sequence_id=";".join(name for name, _ in records),
                   sequence_sha256=digest, length=length, n_contigs=len(records),
                   duplicate_of=hashes.get(digest, ""))
        hashes.setdefault(digest, pid)
    return rows, sequences


def normalize_features(path, index, sequences, assignments):
    if path is None:
        return []
    rows = read_tsv(path, ("plasmid_id", "feature_id", "start", "end", "strand"))
    by_id = {r["plasmid_id"]: r for r in index}
    seen, normalized = set(), []
    for row in rows:
        pid = row["plasmid_id"]
        if pid not in by_id:
            raise ValueError(f"Annotation references unknown plasmid {pid}")
        key = (pid, row["feature_id"])
        if not key[1] or key in seen:
            raise ValueError(f"Duplicate or empty feature ID: {key}")
        seen.add(key)
        source = row.get("source_sequence_id")
        records = dict(sequences[pid])
        if not source and len(records) == 1:
            source = next(iter(records))
        if source not in records:
            raise ValueError(f"{key}: source_sequence_id required for multi-contig plasmids")
        start, end = int(row["start"]), int(row["end"])
        if not 1 <= start <= end <= len(records[source]) or row["strand"] not in {"+", "-", "."}:
            raise ValueError(f"{key}: invalid coordinates (use 1-based inclusive contig coordinates)")
        for field in ("identity", "coverage"):
            if row.get(field) and not 0 <= float(row[field]) <= 100:
                raise ValueError(f"{key}: {field} must be a percentage in [0,100]")
        record = by_id[pid]
        if row.get("isolate_id") and row["isolate_id"] != record["isolate_id"]:
            raise ValueError(f"{key}: annotation isolate mismatch")
        if row.get("sequence_sha256") and row["sequence_sha256"] != record["sequence_sha256"]:
            raise ValueError(f"{key}: annotation checksum mismatch")
        eligible = bool(row.get("amr_gene") and row.get("hit_class", "").lower() in {"strict", "perfect", "curated"}
                        and row.get("annotation_confidence", "").lower() == "high"
                        and row.get("database_name") and row.get("database_version")
                        and row.get("annotation_engine"))
        row.update(isolate_id=record["isolate_id"], plasmid_unit=assignments.get(pid, ""),
                   sequence_sha256=record["sequence_sha256"], source_sequence_id=source,
                   start=start, end=end, headline_eligible=str(eligible).lower())
        if record["quality_status"] != "rejected":
            normalized.append(row)
    return normalized
