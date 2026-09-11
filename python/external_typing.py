"""Optional checksum-linked cross-references to external typing/taxonomy calls.

These fields (e.g. MOB-suite cluster IDs, COPLA plasmid taxonomic units, predicted
host range, PlasmidFinder Inc-types, PLSDB nearest-reference match) are opaque external
identifiers. PlasmiCord does not compute, validate or interpret them, and they never
influence quality/confidence tiers in quality.py. PlasmidFinder and MOB-typer (already
imported) use independently curated replicon databases and may legitimately disagree;
that disagreement is useful review evidence, not something this module reconciles.
"""
import re
from .contracts import checksum, read_tsv

EXTERNAL_TYPING_FIELDS = ("plasmid_id sequence_sha256 evidence_source external_tool external_tool_version "
                          "mob_primary_cluster_id mob_secondary_cluster_id mob_cluster_distance_definition "
                          "ptu_assignment ptu_confidence predicted_host_range_overall_rank "
                          "predicted_host_range_overall_name associated_pmids "
                          "predicted_transmissibility_score predicted_transmissibility_call "
                          "predicted_transmissibility_tool predicted_transmissibility_tool_version "
                          "plasmidfinder_inc_types plasmidfinder_identity "
                          "predicted_classification_score predicted_classification_call "
                          "predicted_classification_tool predicted_classification_tool_version "
                          "plsdb_nearest_accession plsdb_nearest_distance plsdb_nearest_host "
                          "evidence_sha256").split()

TRANSMISSIBILITY_CALLS = {"conjugative", "mobilizable", "non-mobilizable", "uncertain"}
# Matches quality.py's --quality-evidence `classification` vocabulary for consistency; this is a
# different evidentiary axis (identity/classification confidence, e.g. PlasFlow/Plasmer) from the
# transmissibility fields above (transfer-potential confidence, e.g. PlasTrans) -- never conflate them.
CLASSIFICATION_CALLS = {"plasmid", "chromosome", "uncertain"}


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
                or not row["evidence_source"] or not row["external_tool"]):
            raise ValueError(f"External typing must uniquely match a candidate checksum and source, and name an external_tool: {pid}")
        confidence = row.get("ptu_confidence", "")
        if confidence and confidence.lower() not in {"low", "medium", "high"}:
            try:
                value = float(confidence)
            except ValueError:
                raise ValueError(f"{pid}: ptu_confidence must be low/medium/high or a number in [0,100]") from None
            if not 0 <= value <= 100:
                raise ValueError(f"{pid}: ptu_confidence must be low/medium/high or a number in [0,100]")
        pmids = row.get("associated_pmids", "")
        if pmids and not all(re.fullmatch(r"\d+", pmid) for pmid in pmids.split(";")):
            raise ValueError(f"{pid}: associated_pmids must be semicolon-separated numeric PubMed IDs")
        score = row.get("predicted_transmissibility_score", "")
        if score:
            try:
                score_value = float(score)
            except ValueError:
                raise ValueError(f"{pid}: predicted_transmissibility_score must be a number in [0,1]") from None
            if not 0 <= score_value <= 1:
                raise ValueError(f"{pid}: predicted_transmissibility_score must be a number in [0,1]")
        call = row.get("predicted_transmissibility_call", "")
        if call and call.lower() not in TRANSMISSIBILITY_CALLS:
            raise ValueError(f"{pid}: predicted_transmissibility_call must be one of {sorted(TRANSMISSIBILITY_CALLS)}")
        identity = row.get("plasmidfinder_identity", "")
        if identity:
            try:
                identity_value = float(identity)
            except ValueError:
                raise ValueError(f"{pid}: plasmidfinder_identity must be a number in [0,100]") from None
            if not 0 <= identity_value <= 100:
                raise ValueError(f"{pid}: plasmidfinder_identity must be a number in [0,100]")
        classification_score = row.get("predicted_classification_score", "")
        if classification_score:
            try:
                classification_score_value = float(classification_score)
            except ValueError:
                raise ValueError(f"{pid}: predicted_classification_score must be a number in [0,1]") from None
            if not 0 <= classification_score_value <= 1:
                raise ValueError(f"{pid}: predicted_classification_score must be a number in [0,1]")
        classification_call = row.get("predicted_classification_call", "")
        if classification_call and classification_call.lower() not in CLASSIFICATION_CALLS:
            raise ValueError(f"{pid}: predicted_classification_call must be one of {sorted(CLASSIFICATION_CALLS)}")
        plsdb_distance = row.get("plsdb_nearest_distance", "")
        if plsdb_distance:
            try:
                plsdb_distance_value = float(plsdb_distance)
            except ValueError:
                raise ValueError(f"{pid}: plsdb_nearest_distance must be a number in [0,1]") from None
            if not 0 <= plsdb_distance_value <= 1:
                raise ValueError(f"{pid}: plsdb_nearest_distance must be a number in [0,1]")
        row["evidence_sha256"] = evidence_sha256
        typing[pid] = row
    return typing
