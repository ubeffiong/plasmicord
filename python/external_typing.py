"""Optional checksum-linked cross-references to external typing/taxonomy calls.

These fields (e.g. MOB-suite cluster IDs, COPLA plasmid taxonomic units, predicted
host range) are opaque external identifiers. PlasmiCord does not compute, validate
or interpret them, and they never influence quality/confidence tiers in quality.py.
"""
import re
from .contracts import checksum, read_tsv

EXTERNAL_TYPING_FIELDS = ("plasmid_id sequence_sha256 evidence_source external_tool external_tool_version "
                          "mob_primary_cluster_id mob_secondary_cluster_id mob_cluster_distance_definition "
                          "ptu_assignment ptu_confidence predicted_host_range_overall_rank "
                          "predicted_host_range_overall_name associated_pmids evidence_sha256").split()


def load_external_typing(path, index):
    if not path:
        return {}
    rows = read_tsv(path, ("plasmid_id", "sequence_sha256", "evidence_source", "external_tool"))
    by_id = {r["plasmid_id"]: r for r in index}
    evidence_sha256 = checksum(path)
    typing = {}
    for row in rows:
        pid = row["plasmid_id"]
        if (pid not in by_id or pid in typing or row["sequence_sha256"] != by_id[pid]["sequence_sha256"]
                or not row["evidence_source"]):
            raise ValueError(f"External typing must uniquely match a candidate checksum and source: {pid}")
        confidence = row.get("ptu_confidence", "")
        if confidence and confidence.lower() not in {"low", "medium", "high"}:
            try:
                value = float(confidence)
            except ValueError:
                raise ValueError(f"{pid}: ptu_confidence must be low/medium/high or a number in [0,100]")
            if not 0 <= value <= 100:
                raise ValueError(f"{pid}: ptu_confidence must be low/medium/high or a number in [0,100]")
        pmids = row.get("associated_pmids", "")
        if pmids and not all(re.fullmatch(r"\d+", pmid) for pmid in pmids.split(";")):
            raise ValueError(f"{pid}: associated_pmids must be semicolon-separated numeric PubMed IDs")
        row["evidence_sha256"] = evidence_sha256
        typing[pid] = row
    return typing
