"""Native output adapters preserve candidate boundaries and producer provenance."""
import json
import re
from pathlib import Path
from .contracts import fasta_records, read_tsv, write_tsv, checksum, metadata, validate_manifest, resolve_path, INDEX_FIELDS

ADAPTERS = ("generic-fasta", "mob-recon", "flye", "unicycler", "plasbench")


def safe_id(text):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text).strip("._-")


def selected_records(path, selectors):
    records = fasta_records(path)
    if not selectors:
        raise ValueError("Explicit --record-ids (or record_ids in --samples) is required; circular contigs are not automatically plasmids")
    keep = set(selectors.split(","))
    if keep - {name for name, _ in records}:
        raise ValueError(f"Selected record IDs missing from {path}: {sorted(keep - {n for n, _ in records})}")
    return [(name, seq) for name, seq in records if name in keep]


def native_candidates(adapter, source, selectors="", group_contigs=False):
    source = Path(source).resolve()
    if adapter == "mob-recon":
        files = sorted(source.glob("plasmid_*.fasta"))
        report = source / "contig_report.txt"
        if not report.is_file():
            raise ValueError("MOB-recon directory must include contig_report.txt")
        assignments = read_tsv(report, ('contig_id', 'molecule_type'))
        calls = {r['contig_id']: r['molecule_type'].lower() for r in assignments}
        if len(calls) != len(assignments):
            raise ValueError("Duplicate contig assignments in MOB-recon report")
        assigned = set()
        for file in files:
            records = fasta_records(file)
            if assigned & {n for n, seq in records}:
                raise ValueError("A MOB-recon contig occurs in multiple candidate files")
            assigned.update(n for n, seq in records)
            if any(calls.get(n) != 'plasmid' for n, seq in records):
                raise ValueError(f"MOB-recon candidate contigs must have plasmid assignments in contig_report.txt: {file}")
            yield file.stem, records, dict(assembly_method="short_read_reconstruction"), [file, report]
    elif adapter in {"flye", "unicycler"}:
        fasta = source / "assembly.fasta"
        graph = source / ("assembly_graph.gfa" if adapter == "flye" else "assembly.gfa")
        circular, coverage = {}, {}
        evidence = [fasta]
        if adapter == "flye":
            info = source / "assembly_info.txt"
            if not info.is_file():
                raise ValueError("Flye import requires assembly_info.txt")
            with info.open(encoding="utf-8") as handle:
                header = handle.readline().lstrip("#").rstrip().split("\t")
                if 'seq_name' not in header or 'circ.' not in header:
                    raise ValueError("Unsupported Flye assembly_info.txt header")
                for line in handle:
                    r = dict(zip(header, line.rstrip().split("\t")))
                    circular[r['seq_name']] = r.get('circ.') == 'Y'
                    coverage[r['seq_name']] = r.get('cov.', '')
            evidence.append(info)
        else:
            with fasta.open(encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith(">"):
                        circular[line[1:].split()[0]] = 'circular=true' in line.lower()
        if graph.is_file():
            evidence.append(graph)
        for name, seq in selected_records(fasta, selectors):
            if adapter == "flye" and name not in circular:
                raise ValueError(f"Flye assembly_info.txt is missing selected record {name}")
            yield name, [(name, seq)], dict(circularity_status="tool_reported" if circular.get(name) else "unresolved",
                       assembly_method="long_read" if adapter == "flye" else "hybrid",
                       assembly_graph_path=str(graph) if graph.is_file() else "", reported_coverage=coverage.get(name, "")), evidence
    elif adapter == "plasbench":
        reports = sorted(source.glob("*selection_report.json"))
        if not reports:
            raise ValueError("Select the PlasBench selected_candidate(s)/sample directory containing a selection_report.json")
        if len(reports) != 1:
            raise ValueError("Ambiguous PlasBench selection reports; import each sample separately")
        report = json.loads(reports[0].read_text(encoding="utf-8"))
        file = source / "candidate.plasmid.fasta"
        if not file.is_file():
            # Older operational selection exports preserve sample/tool in the filename.
            candidates = sorted(source.glob("*.plasmid.fasta"))
            if len(candidates) != 1:
                raise ValueError("PlasBench selected candidate FASTA missing or ambiguous")
            file = candidates[0]
        records = fasta_records(file)
        if len(records) > 1 and not group_contigs and not selectors:
            raise ValueError("PlasBench multi-record FASTA has ambiguous plasmid boundaries: provide record_ids or --group-contigs for ONE candidate")
        groups = [(file.stem, records)] if group_contigs else [(name, [(name, seq)]) for name, seq in (selected_records(file, selectors) if selectors else records)]
        for name, records in groups:
            yield name, records, dict(source_tool=report.get("selected_tool", "unknown"),
                selection_evidence_path=str(reports[0]), selection_method=report.get("selection_type", "unknown"),
                assembly_method=report.get("analysis_track", "unknown")), [file, reports[0]]
    else:
        records = fasta_records(source)
        if group_contigs:
            yield source.stem, records, {}, [source]
        elif len(records) == 1 and not selectors:
            yield records[0][0], records, {}, [source]
        else:
            for name, seq in selected_records(source, selectors):
                yield name, [(name, seq)], {}, [source]


def import_inputs(args):
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Import output must be a fresh directory")
    if args.group_contigs and args.adapter not in {'generic-fasta', 'plasbench'}:
        raise ValueError('--group-contigs applies only to generic-fasta and plasbench')
    if args.samples:
        samples = read_tsv(args.samples, ("isolate_id", "input_path"))
        for sample in samples:
            sample["input_path"] = str(resolve_path(args.samples.resolve().parent, sample["input_path"]))
    else:
        if not args.input or not args.isolate_id:
            raise ValueError("Provide --samples or both --input and --isolate-id")
        samples = [dict(isolate_id=args.isolate_id, input_path=str(args.input), record_ids=args.record_ids or "")]
    seen, candidates, sources = set(), [], []
    for sample in samples:
        iso = sample["isolate_id"]
        if not iso or iso in seen:
            raise ValueError("Import sample IDs must be unique")
        seen.add(iso)
        for name, records, extra, evidence in native_candidates(args.adapter, sample['input_path'], sample.get('record_ids', ''), args.group_contigs):
            pid = safe_id(f"{iso}__{name}")
            imported = dict(isolate_id=iso, plasmid_id=pid, fasta_path=f"plasmids/{pid}.fasta",
                source_type="plasbench" if args.adapter == "plasbench" else "longread" if args.adapter == "flye" else "hybrid" if args.adapter == "unicycler" else "precomputed",
                source_tool=args.adapter, source_tool_version=sample.get("source_tool_version") or args.tool_version,
                sequencing_technology=sample.get("sequencing_technology", ""))
            imported.update(extra)
            candidates.append((imported, records))
            sources.extend(dict(path=str(p), sha256=checksum(p)) for p in evidence)
    ids = [r["plasmid_id"].casefold() for r, _ in candidates]
    if len(ids) != len(set(ids)):
        raise ValueError("Imported names collide after filename normalization")
    out.mkdir(parents=True, exist_ok=True)
    (out / "plasmids").mkdir()
    for row, records in candidates:
        (out / row["fasta_path"]).write_text("".join(f">{name}\n{seq}\n" for name, seq in records), encoding="utf-8")
    extra_fields = ["reported_coverage", "selection_evidence_path", "selection_method"]
    write_tsv(out / "manifest.tsv", INDEX_FIELDS + extra_fields, [r for r, _ in candidates])
    meta_fields = ["isolate_id", "chromosomal_cluster", "date", "location", "organism"]
    write_tsv(out / "metadata.tsv", meta_fields, samples)
    # Every adapter enters the exact same technical gate before its package is accepted.
    validated, _ = validate_manifest(out / "manifest.tsv", metadata(out / "metadata.tsv"), "precomputed", 1)
    write_tsv(out / "validation.tsv", INDEX_FIELDS, validated)
    (out / "import_provenance.json").write_text(json.dumps(dict(schema_version="1.0", adapter=args.adapter,
        candidate_count=len(candidates), inputs=sources, plasmid_identity="candidate only; evaluate before inference"), indent=2), encoding="utf-8")
    print(f"Imported {len(candidates)} candidates: {out / 'manifest.tsv'}")
