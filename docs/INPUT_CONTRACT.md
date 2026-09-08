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
The last three paths are retained as provenance; they are not automatically analysed.
Use the global `--features` option to import a normalized annotation table.

Circularity vocabulary: `confirmed`, `assembly_supported`, `tool_reported`, `linear`,
`unresolved`. These are submitter-provided observations. The framework never upgrades
tool-reported circularity to confirmation.

Quality vocabulary: `high_confidence`, `moderate_confidence`, `low_confidence`, `uncertain`,
`rejected`. Submitted status is preserved as `declared_quality_status`. The current
technical gate reports `uncertain`, `low_confidence` or `rejected`; it cannot establish
high or moderate biological confidence. Fragmentation, ambiguous bases and undocumented
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
database_name database_version sequence_sha256
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
It is not an orthology or validated core-genome call. Module completeness remains
unresolved. Supplied KO/module/domain fields and provenance are available in gene
details; automated orthology and module analysis are not implemented.

## PlasBench interoperability

Select desired candidate FASTAs and create the same manifest. The `plasbench` mode
records input origin only, without parsing a native export directory. PlasBench's
current protein TSV uses zero-based half-open coordinates; convert `start + 1`, retain
`end`, map `sequence_id` to `source_sequence_id`, and map `gene/product/category`
to `gene_symbol/product_name/functional_category`. Do not blindly copy its annotation
digest: the checksum definitions differ. Project-specific recovery scores remain
upstream evidence and are not treated as transmission confidence.
