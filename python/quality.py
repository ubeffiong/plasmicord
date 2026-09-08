"""Auditable evidence-based quality tiers; rules are research heuristics, not probabilities."""
from collections import Counter
import math
from pathlib import Path
from .contracts import checksum, read_tsv

POLICY_VERSION = "1.0"
QUALITY_FIELDS = "plasmid_id sequence_sha256 quality_status classification replicon_type mobility_class sequence_entropy ambiguous_fraction terminal_overlap_bp graph_closure read_breadth mean_depth chromosome_fraction evidence_source evidence_sha256 policy_version quality_warnings".split()


def entropy(sequence):
    counts = Counter(sequence)
    total = len(sequence)
    return -sum((n / total) * math.log2(n / total) for base, n in counts.items() if base in "ACGT") if total else 0


def terminal_overlap(sequence, minimum=20, maximum=1000):
    for length in range(min(maximum, len(sequence)//2), minimum-1, -1):
        if sequence[:length] == sequence[-length:] and set(sequence[:length]) <= set("ACGT"):
            return length
    return 0


def graph_closure(path, records):
    """Only an exact-sequence single-segment GFA self-link supplies closure support."""
    if not path:
        return "not_evaluated"
    graph = Path(path)
    if not graph.is_file():
        raise ValueError(f"Assembly graph missing: {graph}")
    seqs, closed = {}, set()
    with graph.open(encoding="utf-8") as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if f[0] == 'S' and len(f) >= 3:
                seqs[f[1]] = f[2].upper()
            if f[0] == 'L' and len(f) >= 6 and f[1] == f[3] and f[2] == f[4]:
                closed.add(f[1])
    if len(records) == 1:
        name, seq = records[0]
        if seqs.get(name) == seq and name in closed:
            return "assembly_supported"
    return "unresolved"  # absence of a self-link is not proof of linearity


def load_evidence(path, index):
    if not path:
        return {}
    rows = read_tsv(path, ("plasmid_id", "sequence_sha256", "evidence_source"))
    by_id = {r["plasmid_id"]: r for r in index}
    evidence = {}
    for row in rows:
        pid = row['plasmid_id']
        if pid not in by_id or pid in evidence or row['sequence_sha256'] != by_id[pid]['sequence_sha256'] or not row['evidence_source']:
            raise ValueError(f"Quality evidence must uniquely match a candidate checksum and source: {pid}")
        for key in ('read_breadth', 'chromosome_fraction'):
            if row.get(key) and not 0 <= float(row[key]) <= 1:
                raise ValueError(f"{key} must be a fraction in [0,1]")
        if row.get('mean_depth') and (not math.isfinite(float(row['mean_depth'])) or float(row['mean_depth']) < 0):
            raise ValueError("mean_depth must be finite and nonnegative")
        if row.get('classification', '') not in {'', 'plasmid', 'chromosome', 'uncertain'}:
            raise ValueError("classification must be plasmid, chromosome or uncertain")
        row['evidence_sha256'] = checksum(path)
        evidence[pid] = row
    return evidence


def assess(index, sequences, typing=None, evidence=None):
    output = []
    for row in index:
        pid = row['plasmid_id']
        records = sequences[pid]
        seq = ''.join(seq for _, seq in records)
        obs = (evidence or {}).get(pid, {})
        typ = (typing or {}).get(pid, {})
        replicon = typ.get('rep_type(s)', '')
        marker = any(typ.get(k, '').strip().lower() not in {'', '-', 'na', 'n/a', 'none'}
                     for k in ('rep_type(s)', 'relaxase_type(s)', 'orit_type(s)'))
        classification = obs.get('classification') or ('plasmid_supported' if marker else 'uncertain')
        warnings = [w for w in row.get('quality_warnings', '').split(';') if w and w != 'plasmid_identity_not_independently_validated']
        ent = entropy(seq)
        ambiguous = sum(b not in 'ACGT' for b in seq) / len(seq)
        overlap = terminal_overlap(seq) if len(records) == 1 else 0
        closure = graph_closure(row.get('assembly_graph_path'), records)
        if classification not in {'plasmid', 'plasmid_supported'}:
            warnings.append('plasmid_identity_unresolved')
        if ent < 1.2:
            warnings.append('low_sequence_complexity')
        if overlap:
            warnings.append('duplicated_terminal_overlap_requires_review')
        if closure == 'assembly_supported':
            row['circularity_status'] = 'assembly_supported'
        breadth = float(obs['read_breadth']) if obs.get('read_breadth') else None
        depth = float(obs['mean_depth']) if obs.get('mean_depth') else None
        contamination = float(obs['chromosome_fraction']) if obs.get('chromosome_fraction') else None
        if contamination is not None and contamination > .05:
            warnings.append('chromosomal_contamination_evidence')
        if breadth is not None and breadth < .95:
            warnings.append('incomplete_read_support')
        if depth is not None and depth < 10:
            warnings.append('low_read_depth')
        status = 'uncertain'
        if classification in {'plasmid', 'plasmid_supported'}:
            status = 'moderate_confidence'
        # High tier requires independently supplied identity/read/contamination evidence and sequence-bound graph support.
        if classification == 'plasmid' and closure == 'assembly_supported' and breadth is not None and breadth >= .95 and depth is not None and depth >= 10 and contamination is not None and contamination <= .01:
            status = 'high_confidence'
        if ent < 1.2 or overlap or ambiguous > .01 or row['quality_status'] == 'low_confidence' or any(w in warnings for w in ('chromosomal_contamination_evidence','incomplete_read_support','low_read_depth')):
            status = 'low_confidence'
        if row['quality_status'] == 'rejected' or classification == 'chromosome':
            status = 'rejected'
            if classification == 'chromosome':
                warnings.append('chromosome_classification_evidence')
        row.update(quality_status=status, quality_warnings=';'.join(sorted(set(warnings))))
        output.append(dict(plasmid_id=pid, sequence_sha256=row['sequence_sha256'], quality_status=status,
            classification=classification, replicon_type=replicon, mobility_class=typ.get('predicted_mobility', 'not_evaluated'),
            sequence_entropy=ent, ambiguous_fraction=ambiguous, terminal_overlap_bp=overlap, graph_closure=closure,
            read_breadth=obs.get('read_breadth', ''), mean_depth=obs.get('mean_depth', ''), chromosome_fraction=obs.get('chromosome_fraction', ''),
            evidence_source=obs.get('evidence_source', ''), evidence_sha256=obs.get('evidence_sha256', ''),
            policy_version=POLICY_VERSION, quality_warnings=row['quality_warnings']))
    return output
