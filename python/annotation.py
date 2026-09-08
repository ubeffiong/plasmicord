"""Connect the independent annotation package to PlasmiCord's candidate contract."""
import json
from pathlib import Path
from plasmid_annotation import annotate, to_plasbench
from plasmid_annotation.engine import database_identity
from .contracts import FEATURE_FIELDS, resolve_path, write_tsv


def load_config(path, profile="essential"):
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError(f"{path}: annotation config must be a JSON object")
    if profile == "essential" and any(not config.get(k) for k in ('genes', 'amr', 'mobility')):
        raise ValueError("The essential annotation profile requires genes, amr and mobility configurations; use custom to run a subset")
    for stage in ('genes', 'amr', 'mobility'):
        spec = config.get(stage)
        if spec:
            if not isinstance(spec, dict):
                raise ValueError(f"{stage}: configuration must be a JSON object")
            if not spec.get('database') or not spec.get('database_version'):
                raise ValueError(f"{stage} requires a database and database_version")
            spec['database'] = str(resolve_path(path.parent, spec['database']))
    return config


def annotate_candidates(index, sequences, config_path, out, cache_dir=None, threads=1, profile="essential"):
    config = load_config(config_path, profile)
    if threads < 1:
        raise ValueError("threads must be positive")
    # Fingerprint actual database contents once per run, not once per candidate.
    fingerprints = {stage: database_identity(spec['database'], spec['database_version'])
                    for stage, spec in config.items() if stage in {'genes','amr','mobility'} and spec}
    rows, statuses, typing, pb = [], [], {}, []
    for candidate in index:
        if candidate['quality_status'] == 'rejected':
            continue
        pid = candidate['plasmid_id']
        print(f"[annotation] {pid}", flush=True)
        features, evaluated, profile_row = annotate(sequences[pid], config, Path(out)/'annotation'/pid, cache_dir, threads, fingerprints)
        for f in features:
            f.update(plasmid_id=pid, isolate_id=candidate['isolate_id'], sequence_sha256=candidate['sequence_sha256'])
            rows.append(f)
        statuses.extend(dict(s, plasmid_id=pid) for s in evaluated)
        typing[pid] = profile_row
        pb.extend(dict(f, plasmid_id=pid, isolate_id=candidate['isolate_id']) for f in to_plasbench(features))
    features_file = Path(out)/'automatic_features.tsv'
    write_tsv(features_file, FEATURE_FIELDS, rows)
    (Path(out)/'annotation_provenance.json').write_text(json.dumps(statuses, indent=2), encoding='utf-8')
    write_tsv(Path(out)/'annotation_status.tsv', ['plasmid_id','stage','status','engine','engine_version','database_version','database_sha256','cache','cache_key','n_features'], statuses)
    write_tsv(Path(out)/'plasbench_proteins.tsv', ['plasmid_id','isolate_id','sequence_id','start','end','strand','feature_id','gene','product','category','dbxref','source','version','confidence'], pb)
    (Path(out)/'mobility_typing.json').write_text(json.dumps(typing, indent=2), encoding='utf-8')
    return features_file, statuses, typing
