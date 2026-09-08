"""Versioned local tool execution, normalized features and checksum-verified cache."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.parse import unquote

SCHEMA_VERSION = "1.0"


def digest_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def database_identity(path, version):
    """Hash database contents, not just a mutable path or a user-supplied label."""
    path = Path(path).resolve()
    if not path.is_dir() or not version:
        raise ValueError("Annotation requires a database directory and explicit database_version")
    records = []
    for file in sorted(path.rglob("*")):
        if file.is_file() and '__pycache__' not in file.parts and file.suffix not in {'.pyc', '.log', '.lock'}:
            records.append([str(file.relative_to(path)), digest_file(file)])
    if not records:
        raise ValueError(f"Empty database: {path}")
    return hashlib.sha256(json.dumps(records, separators=(',', ':')).encode()).hexdigest()


def category(gene, product):
    text = (gene + " " + product).lower()
    if any(x in text for x in ("mercur", "arsenic", "arsenate", "copper resistance", "cadmium", "metal resistance")):
        return "metal_resistance"
    if any(x in text for x in ("relaxase", "conjug", "mobilization", "mobilisation")):
        return "mobility"
    if any(x in text for x in ("replication initi", "replication protein")):
        return "replication"
    if any(x in text for x in ("transpos", "integrase", "insertion sequence")):
        return "mobile_element"
    if any(x in text for x in ("partition", "toxin-antitoxin", "maintenance")):
        return "maintenance"
    if "hypothetical" in text:
        return "hypothetical"
    return "other"  # specialized AMR evidence, not a name substring, establishes ARG calls


def parse_gff(path, engine, version, database, database_version):
    rows = []
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("##FASTA"):
                break
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF row in {path}")
            if fields[2] != "CDS":
                continue
            attrs = dict((k, unquote(v)) for item in fields[8].split(";") if "=" in item for k, v in [item.split("=", 1)])
            gene, product = attrs.get("gene", ""), attrs.get("product", "hypothetical protein")
            rows.append(dict(source_sequence_id=fields[0], feature_id=f"{engine}:{attrs.get('ID') or attrs.get('locus_tag') or len(rows)+1}",
                             gene_symbol=gene, product_name=product, feature_type="CDS", start=int(fields[3]), end=int(fields[4]),
                             strand=fields[6], functional_category=category(gene, product), annotation_confidence="predicted",
                             annotation_engine=engine, annotation_engine_version=version, database_name=database,
                             database_version=database_version, dbxref=attrs.get("Dbxref", ""), hit_class="predicted"))
    return rows


def parse_amrfinder(path, version, database_version, min_identity=90.0, min_coverage=90.0):
    rows = []
    with Path(path).open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        aliases = {"Element symbol":"Gene symbol", "Element name":"Sequence name",
                   "Type":"Element type", "Subtype":"Element subtype",
                   "% Identity to reference":"% Identity to reference sequence",
                   "% Coverage of reference":"% Coverage of reference sequence",
                   "Closest reference accession":"Accession of closest sequence"}
        reader.fieldnames = [aliases.get(name, name) for name in (reader.fieldnames or [])]
        required = {"Contig id", "Start", "Stop", "Strand", "Gene symbol", "Sequence name", "Element type", "Method", "% Identity to reference sequence", "% Coverage of reference sequence"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Unsupported AMRFinderPlus output schema; expected nucleotide coordinates and identity/coverage")
        for i, r in enumerate(reader, 1):
            def percentage(value):
                return None if value in {'NA', '', 'N/A'} else float(value)
            identity, coverage = percentage(r["% Identity to reference sequence"]), percentage(r["% Coverage of reference sequence"])
            method, kind, subtype = r["Method"], r["Element type"], r.get("Element subtype", "")
            partial = "PARTIAL" in method or "INTERNAL_STOP" in method or "POINT" in method
            high = (not partial and r.get("Scope", "core").lower() == "core" and
                    identity is not None and coverage is not None and identity >= min_identity and coverage >= min_coverage)
            functional = "amr" if kind == "AMR" else "metal_resistance" if "METAL" in subtype.upper() else "virulence" if kind == "VIRULENCE" else "stress_response"
            start, end = sorted((int(r["Start"]), int(r["Stop"])))
            rows.append(dict(source_sequence_id=r["Contig id"], feature_id=f"amrfinder:{i}", gene_symbol=r["Gene symbol"],
                             product_name=r["Sequence name"], feature_type="AMR_match" if kind == "AMR" else "functional_match",
                             start=start, end=end, strand=r["Strand"], functional_category=functional,
                             functional_subcategory=subtype, amr_gene=r["Gene symbol"] if kind == "AMR" else "",
                             drug_class=r.get("Class", ""), resistance_mechanism="", identity=identity if identity is not None else "",
                             coverage=coverage if coverage is not None else "", annotation_scope=r.get("Scope", "unreported"),
                             hit_class="strict" if high else "partial" if partial else "supporting",
                             annotation_confidence="high" if high else "low", annotation_engine="amrfinder",
                             annotation_engine_version=version, database_name="AMRFinderPlus", database_version=database_version,
                             detection_method=method, reference_accession=r.get("Accession of closest sequence", "")))
    return rows


def to_plasbench(rows):
    """Export the existing PlasBench protein contract; convert 1-based to 0-based."""
    return [dict(sequence_id=r["source_sequence_id"], start=int(r["start"])-1, end=int(r["end"]),
                 strand=r["strand"], feature_id=r["feature_id"], gene=r.get("gene_symbol", ""), product=r.get("product_name", ""),
                 category=r.get("functional_category", "other"), dbxref=r.get("dbxref", ""),
                 source=r["annotation_engine"], version=r.get("annotation_engine_version", ""),
                 confidence=r.get("annotation_confidence", "")) for r in rows if r.get("feature_type") == "CDS"]


def annotate(records, config, workdir, cache_dir=None, threads=1, fingerprints=None):
    """Annotate one candidate. Records are (contig ID, DNA); return features/status/typing.

    Config has named engines: genes (bakta/prokka), amr (amrfinder), mobility (mob_typer).
    Each requires executable, database and database_version. Disabled stages are explicit.
    Caller failures raise; no missing tool silently becomes a negative observation.
    """
    workdir = Path(workdir)
    if threads < 1 or not records or not any(config.get(k) for k in ('genes','amr','mobility')):
        raise ValueError("Annotation requires records, positive threads and at least one configured stage")
    workdir.mkdir(parents=True, exist_ok=True)
    features, statuses, typing = [], [], {}
    canonical = [(f"contig{i+1}", seq) for i, (_, seq) in enumerate(records)]
    mapping = {name: records[i][0] for i, (name, _) in enumerate(canonical)}
    fasta = workdir / "input.fasta"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in canonical), encoding="utf-8")
    sequence_hash = hashlib.sha256("\n".join(seq for _, seq in canonical).encode()).hexdigest()
    for stage in ("genes", "amr", "mobility"):
        spec = config.get(stage)
        if not spec:
            statuses.append(dict(stage=stage, status="not_evaluated", reason="not configured"))
            continue
        engine = spec.get("engine", {"genes": "bakta", "amr": "amrfinder", "mobility": "mob_typer"}[stage])
        allowed = {"genes": {"bakta", "prokka"}, "amr": {"amrfinder"}, "mobility": {"mob_typer"}}
        if engine not in allowed[stage]:
            raise ValueError(f"Unsupported {stage} engine: {engine}")
        executable = shutil.which(spec.get("executable", engine))
        if not executable:
            raise ValueError(f"Required annotation executable unavailable: {engine}")
        v = subprocess.run([executable, "--version"], capture_output=True, text=True, timeout=60, check=True)
        version = (v.stdout or v.stderr).strip()
        if not version:
            raise ValueError(f"{engine} returned no version")
        db = Path(spec["database"]).resolve()
        dbversion = spec["database_version"]
        fingerprint = (fingerprints or {}).get(stage) or database_identity(db, dbversion)
        identity, coverage = float(spec.get("min_identity", 90)), float(spec.get("min_coverage", 90))
        if not 0 <= identity <= 100 or not 0 <= coverage <= 100:
            raise ValueError("AMR thresholds must be percentages within [0,100]")
        provenance = dict(stage=stage, engine=engine, engine_version=version, database_version=dbversion,
                          database_sha256=fingerprint, sequence_sha256=sequence_hash, config=spec, threads=threads,
                          shared_component_version="0.1.0", schema_version=SCHEMA_VERSION)
        key = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
        cached = Path(cache_dir) / f"{key}.json" if cache_dir else None
        result = None
        if cached and cached.is_file():
            document = json.loads(cached.read_text(encoding="utf-8"))
            serialized = json.dumps(document["result"], sort_keys=True)
            if document.get("result_sha256") != hashlib.sha256(serialized.encode()).hexdigest():
                raise ValueError(f"Corrupt annotation cache: {cached}")
            result = document["result"]
        stage_dir = workdir / stage
        stage_dir.mkdir(exist_ok=True)
        if result is None:
            if stage == "genes":
                output_dir = stage_dir / "output"
                if engine == "bakta":
                    command = [executable, "--db", str(db), "--output", str(output_dir), "--prefix", "annotation", "--threads", str(threads), str(fasta)]
                    expected = output_dir / "annotation.gff3"
                else:
                    command = [executable, "--dbdir", str(db), "--outdir", str(output_dir), "--prefix", "annotation", "--cpus", str(threads), str(fasta)]
                    expected = output_dir / "annotation.gff"
            elif stage == "amr":
                expected = stage_dir / "amrfinder.tsv"
                command = [executable, "--nucleotide", str(fasta), "--database", str(db), "--threads", str(threads), "--plus", "--output", str(expected)]
            else:
                expected = stage_dir / "mobtyper.tsv"
                command = [executable, "--infile", str(fasta), "--out_file", str(expected), "--database_directory", str(db), "--num_threads", str(threads)]
            with (stage_dir / "tool.log").open("w", encoding="utf-8") as log:
                subprocess.run(command, stdout=log, stderr=log, check=True, timeout=int(spec.get("timeout_seconds", 7200)))
            if not expected.is_file():
                raise ValueError(f"{engine} completed without expected output: {expected}")
            if stage == "genes":
                result = dict(features=parse_gff(expected, engine, version, engine, dbversion), typing={})
            elif stage == "amr":
                result = dict(features=parse_amrfinder(expected, version, dbversion, identity, coverage), typing={})
            else:
                with expected.open(encoding="utf-8") as handle:
                    parsed = list(csv.DictReader(handle, delimiter="\t"))
                if len(parsed) != 1 or 'rep_type(s)' not in parsed[0]:
                    raise ValueError("Expected a single candidate MOB-typer report")
                result = dict(features=[], typing=parsed[0])
            result.update(raw_output_sha256=digest_file(expected), command=command)
            if cached:
                cached.parent.mkdir(parents=True, exist_ok=True)
                content = dict(result=result, result_sha256=hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest())
                with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=cached.parent, delete=False) as tmp:
                    json.dump(content, tmp)
                    temporary = tmp.name
                os.replace(temporary, cached)
            cache_status = "miss" if cached else "disabled"
        else:
            cache_status = "reused"
        for row in result["features"]:
            if row["source_sequence_id"] not in mapping:
                raise ValueError("Annotation contig does not match input; refusing invalid cached/tool output")
            features.append(dict(row, source_sequence_id=mapping[row["source_sequence_id"]]))
        typing.update(result["typing"])
        status = dict(provenance, status="complete", cache=cache_status, cache_key=key,
                      n_features=len(result["features"]), raw_output_sha256=result["raw_output_sha256"], command=result["command"])
        (stage_dir / "provenance.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
        statuses.append(status)
    return features, statuses, typing
