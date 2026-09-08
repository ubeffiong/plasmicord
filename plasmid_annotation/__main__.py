"""Independent annotation CLI with a PlasBench-compatible protein export."""
import argparse
import csv
import json
from pathlib import Path
from .engine import annotate, to_plasbench


def main():
    parser = argparse.ArgumentParser(description="Shared local plasmid annotation; no PlasmiCord/PlasBench runtime required")
    parser.add_argument('--fasta', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--cache-dir', type=Path)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        parser.error('output must be a fresh directory')
    records = []
    for line in args.fasta.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            records.append([name, ''])
        elif line.strip():
            if not records:
                parser.error('FASTA header required')
            records[-1][1] += line.strip().upper()
    if not records or len({n for n,s in records}) != len(records) or any(not s or set(s)-set('ACGTRYSWKMBDHVN') for n,s in records):
        parser.error('invalid FASTA or duplicate IDs')
    config = json.loads(args.config.read_text(encoding='utf-8'))
    for stage in ('genes','amr','mobility'):
        if config.get(stage):
            db = Path(config[stage]['database'])
            config[stage]['database'] = str(db if db.is_absolute() else args.config.resolve().parent/db)
    features, statuses, typing = annotate(records, config, args.out, args.cache_dir, args.threads)
    (args.out/'features.json').write_text(json.dumps(features,indent=2),encoding='utf-8')
    (args.out/'provenance.json').write_text(json.dumps(statuses,indent=2),encoding='utf-8')
    (args.out/'mobility_typing.json').write_text(json.dumps(typing,indent=2),encoding='utf-8')
    proteins = to_plasbench(features)
    fields = ['sequence_id','start','end','strand','feature_id','gene','product','category','dbxref','source','version','confidence']
    with (args.out/'plasbench_proteins.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(proteins)


if __name__ == '__main__':
    main()
