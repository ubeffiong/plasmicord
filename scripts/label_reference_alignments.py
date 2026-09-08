#!/usr/bin/env python3
"""Independent BLAST alignment comparators for real-reference calibration.

These are sequence-relatedness labels, NOT epidemiological transmission truth.
Pairs are restricted to the same BioProject. Whole BioProjects are allocated to
train/validation before alignment; exact sequences shared across splits are excluded.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from python.contracts import read_tsv, write_tsv, fasta_records, checksum


def union_length(intervals):
    total, right = 0, -1
    for start, end in sorted(intervals):
        if end > right:
            total += end-max(start,right)
            right=end
    return total


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--positive-identity',type=float,default=99.5)
    p.add_argument('--positive-coverage',type=float,default=.9)
    p.add_argument('--negative-coverage',type=float,default=.5)
    args=p.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        p.error('output must be fresh')
    if not 0<args.positive_identity<=100 or not 0<=args.negative_coverage<args.positive_coverage<=1:
        p.error('invalid alignment label thresholds')
    args.out.mkdir(parents=True,exist_ok=True)
    index=read_tsv(args.index,('plasmid_id','isolate_id','fasta_path','sequence_sha256'))
    meta={r['isolate_id']:r for r in read_tsv(args.metadata,('isolate_id','bioproject'))}
    assignments={r['plasmid_id']:('train' if int(hashlib.sha256(meta[r['isolate_id']]['bioproject'].encode()).hexdigest(),16)%2==0 else 'validation') for r in index}
    hash_splits={}
    for r in index:
        hash_splits.setdefault(r['sequence_sha256'],set()).add(assignments[r['plasmid_id']])
    excluded=[r['plasmid_id'] for r in index if len(hash_splits[r['sequence_sha256']])>1]
    retained=[r for r in index if r['plasmid_id'] not in excluded]
    fasta=args.out/'all_plasmids.fasta'
    with fasta.open('w') as handle:
        for r in retained:
            path=Path(r['fasta_path']);path=path if path.is_absolute() else args.index.resolve().parent/path
            records=fasta_records(path)
            if len(records)!=1:
                p.error('Reference comparator requires complete single-contig plasmids')
            handle.write(f">{r['plasmid_id']}\n{records[0][1]}\n")
    version=subprocess.check_output(['blastn','-version'],text=True).strip()
    raw=args.out/'blastn.tsv'
    command=['blastn','-query',str(fasta),'-subject',str(fasta),'-task','megablast','-outfmt','6 qseqid sseqid length nident qstart qend sstart send qlen slen',
             '-max_target_seqs','1000','-evalue','1e-10','-num_threads',str(args.threads),'-out',str(raw)]
    subprocess.run(command,check=True)
    hits={}
    for line in raw.read_text().splitlines():
        a,b,length,nident,qs,qe,ss,se,qlen,slen=line.split('\t')
        if a!=b:
            hits.setdefault((a,b),[]).append((int(length),int(nident),min(int(qs),int(qe))-1,max(int(qs),int(qe)),min(int(ss),int(se))-1,max(int(ss),int(se)),int(qlen),int(slen)))
    rows=[]
    raw_sha256=checksum(raw)
    for a,b in itertools.combinations(retained,2):
        group=meta[a['isolate_id']]['bioproject']
        if not group or group!=meta[b['isolate_id']]['bioproject'] or a['isolate_id']==b['isolate_id']:
            continue
        forward=hits.get((a['plasmid_id'],b['plasmid_id']),[])
        reverse=hits.get((b['plasmid_id'],a['plasmid_id']),[])
        def score(hs):
            if not hs:return 0,0
            ident=100*sum(h[1] for h in hs)/sum(h[0] for h in hs)
            cov=min(union_length([(h[2],h[3]) for h in hs])/hs[0][6],union_length([(h[4],h[5]) for h in hs])/hs[0][7])
            return ident,cov
        i1,c1=score(forward);i2,c2=score(reverse)
        identity,coverage=min(i1,i2),min(c1,c2)
        label='positive' if identity>=args.positive_identity and coverage>=args.positive_coverage else 'negative' if max(c1,c2)<args.negative_coverage else 'uncertain'
        rows.append(dict(plasmid_a=a['plasmid_id'],plasmid_b=b['plasmid_id'],label=label,split=assignments[a['plasmid_id']],group_id=group,
                         evidence_source=f"BLASTn:{raw_sha256};NCBI_BioProject:{group}",evidence_type='independent_alignment',
                         reciprocal_identity=identity,reciprocal_coverage=coverage))
    write_tsv(args.out/'labels.tsv',['plasmid_a','plasmid_b','label','split','group_id','evidence_source','evidence_type','reciprocal_identity','reciprocal_coverage'],rows)
    (args.out/'label_provenance.json').write_text(json.dumps(dict(version=version,command=command,excluded_cross_split_duplicates=excluded,
         positive_identity=args.positive_identity,positive_coverage=args.positive_coverage,negative_coverage=args.negative_coverage,
         index_sha256=checksum(args.index),metadata_sha256=checksum(args.metadata),raw_sha256=checksum(raw),
         meaning='Independent sequence alignment comparator, not reviewed transmission truth'),indent=2))
    print(f'Wrote {len(rows)} alignment comparator pairs to {args.out}')


if __name__=='__main__':
    main()
