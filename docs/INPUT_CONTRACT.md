# Standalone input contract, v1

The core begins with reconstructed plasmid candidates, not raw sequencing reads.
No producer-specific code or PlasBench runtime is needed.

## Candidate manifest

Required: `isolate_id`, `plasmid_id`, `fasta_path`. Candidate IDs are unique even
ignoring case, portable filename tokens (letters, numbers, `_`, `-`, `.`),
and cannot use Windows reserved filenames. Metadata IDs must be unique and each
candidate must map to a metadata isolate. Relative paths resolve from the manifest.

Optional provenance: `source_type`, `source_tool`, `source_tool_version`,
`sequencing_technology`, `assembly_method`, `polishing_method`, `circularity_status`,
`quality_status`, `sequence_sha256`, `annotation_path`, `reads_path`, `assembly_graph_path`.
An assembly graph can provide sequence-bound closure evidence. Annotation/read paths
remain provenance; use `--features`, `--annotation-config` or `--quality-evidence` to
supply the corresponding analysis. See [biological quality](QUALITY.md).

Circularity vocabulary: `confirmed`, `assembly_supported`, `tool_reported`, `linear`,
`unresolved`. These are submitter-provided observations. The framework never upgrades
tool-reported circularity to confirmation.

Quality vocabulary: `high_confidence`, `moderate_confidence`, `low_confidence`, `uncertain`,
`rejected`. Submitted status is preserved as `declared_quality_status`. The current
technical gate reports `uncertain`, `low_confidence` or `rejected`; the subsequent
[biological evidence policy](QUALITY.md) may assign moderate or high tiers when its
explicit requirements are met. Fragmentation, ambiguous bases and undocumented
ONT polishing produce warnings. Candidates below `--min-length` (default 200 bases)
or explicitly rejected by the submitter are excluded. All others remain visible with
their limitations; a missing replicon never causes rejection.

Valid FASTA, nonempty contigs, unique identifiers within each candidate, an IUPAC DNA
alphabet and provenance mapping are mandatory. Multi-contig candidates remain separate
records, and distance computation never introduces k-mers across contig joins.
Sequence SHA-256 hashes uppercase contig sequences joined by a newline, without a final
newline, in input order; sequence names and line wrapping are excluded. Duplicate
hashes are flagged but retained across isolates. This does not normalize circular
rotation, reverse complement or reordered contigs.

Malformed technical inputs abort before analysis. Rejected candidates are listed in
`validation.tsv`. Runs require a fresh output folder to prevent stale-result mixing.

## Isolate metadata

Required: `isolate_id`. Optional: `chromosomal_cluster`, `date` (`YYYY-MM-DD`),
`location`, `organism`. Every metadata isolate appears in the network even if it has
no accepted indexed plasmid. This is not a claim that it is plasmid-free.
Unknown, NA, N/A, NONE and blank chromosomal clusters are excluded from discordance
calls. Missing dates/locations remain missing. No link is temporally oriented.
Chromosomal labels must be comparable and namespaced across species or typing schemes.

## Functional feature TSV

Required: `plasmid_id`, `feature_id`, `start`, `end`, `strand`.
Coordinates are **1-based inclusive on the original contig**, not a concatenated
plasmid. `source_sequence_id` is required for multi-contig candidates and must match
a FASTA record. Strand is `+`, `-`, or `.`. Coordinates outside the record are rejected.

Additional normalized columns:

```text
isolate_id plasmid_unit gene_symbol product_name feature_type
functional_category functional_subcategory amr_gene drug_class resistance_mechanism
replicon_type mobility_function ko_id kegg_module cog_category go_terms ec_number
pfam_ids identity coverage hit_class annotation_confidence annotation_engine
database_name database_version sequence_sha256 annotation_engine_version
detection_method reference_accession dbxref
```

Identity/coverage, when provided, are percentages from 0 to 100. Isolate, unit and
sequence digest are derived from validated inputs; supplied isolate/digest mismatches
are errors. Feature IDs are unique within each plasmid. Categories are caller-supplied;
the importer does not infer resistance from a protein name such as `merA`.

Headline ARG eligibility requires `amr_gene`, `hit_class` of strict/perfect/curated,
`annotation_confidence=high`, and nonempty engine/database/version fields. This
filters imported evidence; it does not rerun or validate database thresholds.
Loose, low-confidence and incomplete-provenance calls remain inspectable without
contributing to headline ARG filters. The importer does not combine independent
databases or resolve contradictory alleles; normalize one primary caller upstream.

Unit prevalence = unique candidate plasmids carrying an annotation label divided
by all accepted candidate plasmids in that unit. Counts are deduplicated across
feature copies. `core_or_accessory` is deliberately `observed_in_all` or
`observed_in_subset`, since missing annotation cannot demonstrate biological absence.
Automatic runs also report `n_evaluated_plasmids` and `annotation_coverage` for the
relevant caller stage; imported tables cannot imply completed zero-hit searches.
It is not an orthology or validated core-genome call. Module completeness remains
unresolved. Supplied KO/module/domain fields and provenance are available in gene
details; automated orthology and module analysis are not implemented.

## External typing/taxonomy cross-reference

`--external-typing` accepts an optional checksum-linked TSV carrying identifiers
PlasmiCord does not compute itself: `mob_primary_cluster_id`, `mob_secondary_cluster_id`,
`mob_cluster_distance_definition`, `ptu_assignment`, `ptu_confidence`,
`predicted_host_range_overall_rank`, `predicted_host_range_overall_name`,
`associated_pmids`. Required columns are `plasmid_id`, `sequence_sha256`,
`evidence_source`, `external_tool`. `sequence_sha256` must match the candidate's checksum;
`plasmid_id` must be unique in the file. `ptu_confidence` is `low`/`medium`/`high` or a
number in `[0,100]`. `associated_pmids` is semicolon-separated numeric PubMed IDs.

**These fields are opaque external identifiers.** PlasmiCord does not validate, recompute
or interpret them, and they never influence quality/confidence tiers -- see
[biological quality](QUALITY.md). Output: `typing_crossreference.tsv`.

## Cluster-assignment margin and containment heuristic

`network.edge_evidence.tsv` includes `threshold_margin` (`threshold - minimum_distance`)
per edge, and `plasmid_clusters.tsv` includes `unit_max_internal_distance` and
`unit_threshold_margin` per plasmid unit (blank for singletons). Under complete linkage
the margin is always non-negative; under single linkage it can go negative, exposing a
unit whose members are only connected through chained membership rather than mutual
similarity within the threshold.

`containment_candidates.tsv` flags near-equal-length, high-similarity plasmid pairs using
`--containment-min-ratio`/`--containment-max-ratio`/`--containment-max-distance` (defaults
0.5, 0.95, and `--threshold`). This is a length/similarity heuristic only, **not**
alignment-confirmed containment, and cannot detect size-disparate containment (a small
plasmid nested in a much larger one) -- see [REVIEW.md](REVIEW.md) for why.

## PlasBench interoperability

Use `plasmicord import plasbench` for native selected-candidate exports, or create
the same manifest directly. The `--mode plasbench` option records the input origin.
See [native adapters](IMPORTS.md) and [shared annotation](../shared_annotation/README.md). PlasBench's
current protein TSV uses zero-based half-open coordinates; convert `start + 1`, retain
`end`, map `sequence_id` to `source_sequence_id`, and map `gene/product/category`
to `gene_symbol/product_name/functional_category`. Do not blindly copy its annotation
digest: the checksum definitions differ. Project-specific recovery scores remain
upstream evidence and are not treated as transmission confidence.
