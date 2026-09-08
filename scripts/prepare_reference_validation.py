#!/usr/bin/env python3
"""Copy accession-backed NCBI plasmids into a standalone, auditable validation cohort.

Accepts a generic cohort TSV and local reference layout (sample/reference.fna plus
NCBI sequence_report.jsonl); no producer package is imported. No epidemiological
cluster, date or location is invented from reference sequence similarity.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from python.contracts import fasta_records, write_tsv, checksum


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', type=Path, required=True)
    p.add_argument('--references', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args=p.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        p.error('output must be fresh')
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'plasmids').mkdir()
    cohort=list(csv.DictReader(args.cohort.open(),delimiter='\t'))
    manifest, meta, evidence, sources, excluded=[],[],[],[],[]
    for sample in cohort:
        sid=sample['sample_id']
        folder=args.references/sid
        reference, report=folder/'reference.fna', folder/'sequence_report.jsonl'
        if not reference.is_file() or not report.is_file():
            excluded.append(dict(sample_id=sid,reason='local reference FASTA or NCBI sequence report unavailable'))
            continue
        ncbi=[json.loads(line) for line in report.read_text().splitlines() if line.strip()]
        plasmids={r.get('refseqAccession') or r.get('genbankAccession'):r for r in ncbi if r.get('assignedMoleculeLocationType')=='Plasmid'}
        for name, seq in fasta_records(reference):
            if name not in plasmids:
                continue
            annotation=plasmids[name]
            if annotation.get('assemblyAccession')!=sample['assembly_accession'] or int(annotation['length'])!=len(seq):
                raise ValueError('NCBI sequence report does not match accession/sequence length')
            pid=f'{sid}__{name}'
            path=f'plasmids/{pid}.fasta'
            (args.out/path).write_text(f'>{name}\n{seq}\n')
            digest=hashlib.sha256(seq.encode()).hexdigest()
            manifest.append(dict(isolate_id=sid,plasmid_id=pid,fasta_path=path,source_type='public_reference',source_tool='NCBI_Datasets',
                                 source_tool_version='recorded_in_source_package',assembly_method=sample.get('truth_technology',''),
                                 circularity_status='unresolved',sequence_sha256=digest))
            evidence.append(dict(plasmid_id=pid,sequence_sha256=digest,classification='plasmid',
                                 evidence_source=f"https://www.ncbi.nlm.nih.gov/nuccore/{name}"))
        meta.append(dict(isolate_id=sid,organism=sample.get('organism',''),bioproject=sample.get('bioproject',''),
                         source_study=sample.get('source_study',''),assembly_accession=sample['assembly_accession']))
        sources.append(dict(sample_id=sid,assembly_accession=sample['assembly_accession'],bioproject=sample.get('bioproject',''),
                            reference_sha256=checksum(reference),sequence_report_sha256=checksum(report)))
    if not manifest:
        p.error('No NCBI-designated plasmid sequences found in the supplied cohort')
    write_tsv(args.out/'manifest.tsv',list(manifest[0]),manifest)
    write_tsv(args.out/'metadata.tsv',['isolate_id','organism','bioproject','source_study','assembly_accession'],meta)
    write_tsv(args.out/'quality_evidence.tsv',['plasmid_id','sequence_sha256','classification','evidence_source'],evidence)
    (args.out/'source_provenance.json').write_text(json.dumps(dict(kind='real_public_reference_sequences',sources=sources,
        n_cohort_samples=len(cohort),excluded_samples=excluded,
        cohort_sha256=checksum(args.cohort),limitation='NCBI plasmid classification is reference evidence; no transmission labels or verified closure are supplied'),indent=2))
    print(f'{len(manifest)} public plasmids from {len(meta)} isolates copied to {args.out}')


if __name__=='__main__':
    main()
