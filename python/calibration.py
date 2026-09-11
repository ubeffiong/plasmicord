"""Calibrate a PU distance threshold against independent, explicitly scoped labels."""
import html
import json
import math
from pathlib import Path
from .cluster_plasmids import read_matrix, cluster_single, cluster_complete
from .contracts import read_tsv, write_tsv, checksum


def metrics(rows, predicted):
    tp = sum(r['label'] == 'positive' and predicted[i] for i, r in enumerate(rows))
    fp = sum(r['label'] == 'negative' and predicted[i] for i, r in enumerate(rows))
    fn = sum(r['label'] == 'positive' and not predicted[i] for i, r in enumerate(rows))
    tn = sum(r['label'] == 'negative' and not predicted[i] for i, r in enumerate(rows))
    precision = tp/(tp+fp) if tp+fp else None
    recall = tp/(tp+fn) if tp+fn else None
    specificity = tn/(tn+fp) if tn+fp else None
    return dict(tp=tp, fp=fp, fn=fn, tn=tn, precision=precision, recall=recall, specificity=specificity,
                f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0,
                balanced_accuracy=(recall+specificity)/2 if recall is not None and specificity is not None else None)


def calibrate(args):
    ids, positions, matrix = read_matrix(args.matrix)
    index = read_tsv(args.index, ('plasmid_id', 'isolate_id', 'sequence_sha256'))
    by_id = {r['plasmid_id']: r for r in index}
    if len(by_id) != len(index) or set(ids) != set(by_id):
        raise ValueError("Calibration index and matrix must contain the same unique candidates")
    if any(len(row['sequence_sha256']) != 64 or any(c not in '0123456789abcdef' for c in row['sequence_sha256']) for row in index):
        raise ValueError("Calibration index requires lowercase SHA-256 sequence digests")
    run_path = args.matrix.resolve().parent / 'run_provenance.json'
    run_provenance = {}
    if run_path.is_file():
        run_provenance = json.loads(run_path.read_text(encoding='utf-8'))
        if run_provenance.get('status') != 'complete' or run_provenance.get('distance_engine') != args.engine:
            raise ValueError("Calibration engine must match a completed source run")
        expected = run_provenance.get('output_sha256', {}).get(args.matrix.name)
        if expected and expected != checksum(args.matrix):
            raise ValueError("Calibration matrix checksum differs from source run")
        expected_index = run_provenance.get('output_sha256', {}).get(args.index.name)
        if expected_index and expected_index != checksum(args.index):
            raise ValueError("Calibration index checksum differs from source run")
    rows = read_tsv(args.labels, ('plasmid_a', 'plasmid_b', 'label', 'split', 'group_id', 'evidence_source', 'evidence_type'))
    seen, split_ids, group_splits = set(), {'train':set(), 'validation':set()}, {}
    retained = []
    for row in rows:
        a, b = row['plasmid_a'], row['plasmid_b']
        pair = tuple(sorted((a,b)))
        if a == b or a not in positions or b not in positions or pair in seen:
            raise ValueError("Calibration pairs must be distinct, unique and present in the matrix")
        seen.add(pair)
        if row['label'] not in {'positive','negative','uncertain'} or row['split'] not in split_ids:
            raise ValueError("Use positive/negative/uncertain labels and train/validation splits")
        if not row['evidence_source'] or not row['group_id']:
            raise ValueError("Calibration requires evidence sources and independent grouping IDs")
        if row['evidence_type'] not in {'reviewed_plasmid_sharing', 'experimental', 'independent_alignment', 'synthetic_truth'}:
            raise ValueError("Unsupported label evidence: labels derived from Mash/PU membership cannot calibrate Mash")
        if args.cohort_kind == 'real' and row['evidence_type'] == 'synthetic_truth':
            raise ValueError("Synthetic labels cannot be presented as real-cohort calibration")
        group = row['group_id']
        if group in group_splits and group_splits[group] != row['split']:
            raise ValueError("A group spans training and validation; split by independent study/patient groups")
        group_splits[group] = row['split']
        split_ids[row['split']].update((a,b))
        if row['label'] != 'uncertain':
            retained.append(row)
    for field in ('plasmid_id', 'isolate_id', 'sequence_sha256'):
        left = {by_id[pid][field] for pid in split_ids['train']}
        right = {by_id[pid][field] for pid in split_ids['validation']}
        if '' in left or '' in right or left & right:
            raise ValueError(f"Training/validation leakage or missing values in {field}")
    try:
        thresholds = sorted(set(float(t) for t in args.thresholds.split(',')))
    except ValueError:
        raise ValueError("Thresholds must be a comma-separated list of numbers") from None
    if not thresholds or any(not math.isfinite(t) or not 0 <= t <= 1 for t in thresholds):
        raise ValueError("Thresholds must be finite distances in [0,1]")
    if args.min_groups < 1:
        raise ValueError("min-groups must be positive")
    partitions = {s:[r for r in retained if r['split'] == s] for s in split_ids}
    def predict(split, threshold):
        names = sorted(split_ids[split])
        distances = [[matrix[positions[a]][positions[b]] for b in names] for a in names]
        method = cluster_complete if args.linkage == 'complete' else cluster_single
        unit = {names[i]: c for c, group in enumerate(method(names, distances, threshold)) for i in group}
        return [unit[r['plasmid_a']] == unit[r['plasmid_b']] for r in partitions[split]]
    scores = [dict(split='train', threshold=t, **metrics(partitions['train'], predict('train', t))) for t in thresholds]
    enough = all({r['label'] for r in partitions[s]} == {'positive','negative'} and
                 len({r['group_id'] for r in partitions[s]}) >= args.min_groups for s in partitions)
    chosen = None
    if enough:
        # Holdout metrics never participate in selection.
        chosen = max(scores, key=lambda r:(r['balanced_accuracy'], r['f1'], -r['threshold']))['threshold']
        scores.append(dict(split='validation', threshold=chosen, **metrics(partitions['validation'], predict('validation', chosen))))
    out = args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Calibration output must be a fresh directory")
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out/'calibration_metrics.tsv', ['split','threshold','tp','fp','fn','tn','precision','recall','specificity','f1','balanced_accuracy'], scores)
    write_tsv(out/'calibration_labels.tsv', list(rows[0]) if rows else ['plasmid_a','plasmid_b','label','split','group_id','evidence_source','evidence_type'], rows)
    evidence_types = sorted({r['evidence_type'] for r in retained})
    record = dict(schema_version='1.0', status='evaluated' if enough else 'insufficient_independent_labels',
                  cohort_kind=args.cohort_kind, selected_threshold=chosen, distance_engine=args.engine, linkage=args.linkage,
                  target='plasmid_sequence_relatedness' if 'independent_alignment' in evidence_types else 'reviewed_plasmid_sharing',
                  min_groups_per_split=args.min_groups, label_evidence_types=evidence_types,
                  independent_groups={s:len({r['group_id'] for r in partitions[s]}) for s in partitions},
                  n_pairs={s:len(partitions[s]) for s in partitions}, n_uncertain=sum(r['label']=='uncertain' for r in rows),
                  matrix_sha256=checksum(args.matrix), index_sha256=checksum(args.index), labels_sha256=checksum(args.labels),
                  source_run_sha256=checksum(run_path) if run_path.is_file() else None,
                  source_parameters={k:run_provenance[k] for k in ('k','sketch_size','mash_version') if k in run_provenance},
                  engine_provenance='source_run_verified' if run_provenance else 'submitter_declared',
                  selection_rule='maximize training balanced accuracy, then F1, then smallest threshold',
                  limitation='Sequence similarity and reference alignment agreement do not validate direct transmission. Pair counts are correlated within groups; no population confidence interval is claimed.')
    (out/'calibration.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    narrative = f"# PlasmiCord threshold calibration\n\nStatus: **{record['status']}**\n\nCohort: {args.cohort_kind}; target: {record['target']}.\n\nSelected threshold: {chosen if chosen is not None else 'not selected — insufficient independent positive/negative labels in both splits'}.\n\n{record['limitation']}\n\nTraining and validation do not share candidate IDs, isolate IDs, sequence checksums or groups.\n"
    (out/'CALIBRATION.md').write_text(narrative, encoding='utf-8')
    (out/'CALIBRATION.html').write_text('<!doctype html><meta charset="utf-8"><title>PlasmiCord calibration</title><style>body{font:16px/1.6 system-ui;max-width:1000px;margin:40px auto;padding:20px}pre{white-space:pre-wrap}</style><h1>PlasmiCord threshold calibration</h1><pre>'+html.escape(narrative)+'\n'+html.escape(json.dumps(record,indent=2))+'</pre><p><a href="calibration_metrics.tsv">Metrics TSV</a> · <a href="calibration.json">Provenance JSON</a></p>', encoding='utf-8')
    fields = ('split','threshold','tp','fp','fn','tn','precision','recall','balanced_accuracy')
    table = '<h2>Training sweep and selected holdout result</h2><table cellpadding="8"><tr>'+''.join('<th>'+k+'</th>' for k in fields)+'</tr>'
    for score in scores:
        table += '<tr>'+''.join('<td>'+html.escape(str(score[k]) if score[k] is not None else 'not estimable')+'</td>' for k in fields)+'</tr>'
    with (out/'CALIBRATION.html').open('a', encoding='utf-8') as handle:
        handle.write(table+'</table>')
    print(f"Calibration {record['status']}: {out/'CALIBRATION.html'}")


def add_parser(subs):
    p = subs.add_parser('calibrate', help='Fit a threshold on training labels and evaluate an independent holdout')
    p.add_argument('--matrix', type=Path, required=True)
    p.add_argument('--index', type=Path, required=True)
    p.add_argument('--labels', type=Path, required=True)
    p.add_argument('--cohort-kind', choices=('real','synthetic'), required=True)
    p.add_argument('--engine', choices=('mash','kmer'), required=True)
    p.add_argument('--thresholds', default='0,0.001,0.002,0.005,0.01,0.02,0.05')
    p.add_argument('--linkage', choices=('single','complete'), default='complete')
    p.add_argument('--min-groups', type=int, default=3)
    p.add_argument('--out', type=Path, required=True)
    p.set_defaults(func=calibrate)
