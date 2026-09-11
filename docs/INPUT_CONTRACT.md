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
`associated_pmids`, `predicted_transmissibility_score`, `predicted_transmissibility_call`,
`predicted_transmissibility_tool`, `predicted_transmissibility_tool_version`,
`plasmidfinder_inc_types`, `plasmidfinder_identity`, `predicted_classification_score`,
`predicted_classification_call`, `predicted_classification_tool`,
`predicted_classification_tool_version`, `plsdb_nearest_accession`, `plsdb_nearest_distance`,
`plsdb_nearest_host`. Required columns are `plasmid_id`,
`sequence_sha256`, `evidence_source`, `external_tool`. `sequence_sha256` must match the
candidate's checksum; `plasmid_id` must be unique in the file. `ptu_confidence` is
`low`/`medium`/`high` or a number in `[0,100]`. `associated_pmids` is semicolon-separated
numeric PubMed IDs. `predicted_transmissibility_score` is a number in `[0,1]`;
`predicted_transmissibility_call` is `conjugative`/`mobilizable`/`non-mobilizable`/`uncertain`
-- generic fields for any sequence-based transmissibility classifier's output (e.g.
PlasTrans's codon-usage CNN score), not tied to
one specific tool. `plasmidfinder_inc_types` (semicolon-separated) and
`plasmidfinder_identity` (a number in `[0,100]`) carry
PlasmidFinder's independently-sourced
Inc-type replicon calls. PlasmidFinder and MOB-typer's `rep_type(s)` (already imported) use
different curated databases and may legitimately disagree on the same plasmid; that
disagreement is review evidence, not something PlasmiCord reconciles.
`predicted_classification_score` (a number in `[0,1]`) and `predicted_classification_call`
(`plasmid`/`chromosome`/`uncertain`, matching `--quality-evidence`'s own vocabulary) carry a
composition-based classifier's identity confidence -- e.g.
PlasFlow's or
Plasmer's per-contig probability. This is a **different
evidentiary axis from `predicted_transmissibility_score`**: classification confidence (is this
sequence a plasmid at all) versus transfer-potential confidence (could this plasmid conjugate) --
the two must never be conflated. As with every field in this section, a numeric score here does
not feed `--quality-evidence` or influence PlasmiCord's own quality tiers; only that same tool's
*categorical* call, if separately supplied via `--quality-evidence`'s `classification` field
(see [biological quality](QUALITY.md)), does that.
`plsdb_nearest_accession`, `plsdb_nearest_distance` (a number in `[0,1]`, a Mash-style distance
not a percentage) and `plsdb_nearest_host` carry a candidate's nearest match against
PLSDB (~72,360 curated, dereplicated complete plasmid
sequences derived from NCBI, CC-BY licensed). PlasmiCord does not maintain or bundle a reference
plasmid database itself; a user downloads PLSDB's own Mash-sketch flat file, runs `mash dist`
locally against their accepted candidates, and supplies the result here -- no new PlasmiCord
dependency, same pattern as every other field in this section. Attribute PLSDB per its CC-BY
terms when citing this evidence.

**These fields are opaque external identifiers.** PlasmiCord does not validate, recompute
or interpret them, and they never influence quality/confidence tiers -- see
[biological quality](QUALITY.md). Output: `typing_crossreference.tsv`, written only when
`--external-typing` is supplied, scoped to accepted (non-rejected) candidates.

## Cluster-assignment margin and containment heuristic

`network.edge_evidence.tsv` includes `threshold_margin` (`threshold - minimum_distance`)
per edge, and `plasmid_clusters.tsv` includes `unit_max_internal_distance` and
`unit_threshold_margin` per plasmid unit (blank for singletons). Under complete linkage
the margin is always non-negative; under single linkage it can go negative, exposing a
unit whose members are only connected through chained membership rather than mutual
similarity within the threshold.

`containment_candidates.tsv` is always computed (no flag required) and flags
near-equal-length, high-similarity plasmid pairs using
`--containment-min-ratio`/`--containment-max-ratio`/`--containment-max-distance` (defaults
0.5, 0.95, and `--threshold`; bounds are validated upfront, before any expensive work). This
is a length/similarity heuristic only, **not** alignment-confirmed containment, and cannot
detect size-disparate containment (a small plasmid nested in a much larger one) -- see
[REVIEW.md](REVIEW.md) for why. The default `--containment-max-ratio 0.95` deliberately
excludes equal-length pairs (identity/PU territory, not containment); widening it to `1.0`
would include them.

## Size-corrected distance threshold

`--size-correction-per-percent` (default `0.0`, disabled) optionally loosens the clustering
distance threshold as a plasmid pair's length difference grows, following Scherff et al.
(*Real-time Plasmid Transmission Detection Pipeline*, Microbiol Spectrum 2024,
doi:10.1128/spectrum.02100-24), whose worked example (a 17%-larger plasmid, raw distance 0.003,
only merging once size-corrected) confirms the effective threshold **loosens**, not tightens:

```text
effective_threshold = threshold + size_correction_per_percent * min(size_diff_pct, size_correction_cap_pct)
size_diff_pct = 100 * |len_a - len_b| / max(len_a, len_b)
```

`--size-correction-cap-pct` (default `40.0`, matching the paper) bounds how far the correction
grows. **The paper does not pin down its exact percent-difference denominator convention**;
`max(len_a, len_b)` here is the standard definition, documented as PlasmiCord's own choice, not
asserted as the paper's literal arithmetic. Disabled by default (`0.0`) reproduces clustering
exactly as before -- this is the one option in PlasmiCord that changes clustering results, so it
is opt-in, unlike purely additive diagnostics elsewhere. Affects `plasmid_clusters.tsv`,
`threshold_sensitivity.tsv`, and `network.edge_evidence.tsv` (whose `interpretation` notes when
a link's direct support depends on the correction). `plasmid_clusters.tsv`'s
`unit_threshold_margin` is always reported against the base `--threshold`, never the corrected
one -- a single per-unit margin cannot represent per-pair effective thresholds that vary by
member length.

## Multilayer network export

`--multilayer-network` (opt-in; no files are written when omitted) additionally writes
`network.multilayer.graphml` and `network.multilayer_edges.tsv`. Unlike the default
`network.graphml`/`network.edges.tsv` (which aggregate all shared plasmid units into one
edge per isolate pair), the multilayer export writes **one edge per `(isolate pair,
plasmid unit)`** -- each plasmid unit is naturally one layer -- and tags every edge with
`cluster_relation`: `same_cluster`, `cross_cluster`, or `unknown_cluster` (either isolate's
`chromosomal_cluster` is blank), directly surfacing the relation the discordance feature
already cares about. This is a lightweight, stdlib-only relabeling of edges PlasmiCord
already computes by fields it already has -- not a true node-multiplex/projection framework
(e.g. `pymnet`/`multinet`); see [REVIEW.md](REVIEW.md) for the scoping note.

## PlasAnn interoperability

PlasAnn (Prodigal+BLAST+Infernal) overlaps mostly
with annotation PlasmiCord already gets from wrapped Prokka/Bakta/AMRFinderPlus/MOB-typer,
except for oriT/oriV and transposon calls, which are not otherwise covered. Its output maps
onto the functional feature TSV's existing generic columns without any schema change:
gene coordinates/strand to `start`/`end`/`strand`; oriT/oriV presence to
`functional_category=origin_of_transfer`, `gene_symbol`/`product_name` describing the
element; transposon calls to `functional_category=transposon`; set
`detection_method`/`annotation_engine=PlasAnn` and `reference_accession` to the matched
database entry. Use this mapping only for oriT/oriV/transposon evidence -- PlasAnn is not a
replacement for the wrapped annotation pipeline.

## ISfinder-sequences interoperability

PlasmiCord has no insertion-sequence (IS) detection capability today. A user who BLASTs their
candidates against ISfinder-sequences
(`IS.fna`/`IS.faa`) externally can map hits into `functional_features.tsv` with no schema
change: the IS family/group name to `gene_symbol`, the ISfinder Accession Number to
`reference_accession`, `database_name=ISfinder`, `detection_method`/`annotation_engine` set to
the BLAST tool/version used, and a new documented `functional_category` value:
`insertion_sequence`. **Caveat**: this mirror carries no license and has been stale since 2020
(four commits total) -- treat it as a convenience source only. Cite and verify against the
[ISfinder](https://isfinder.biotoul.fr/) database directly rather than treating this mirror as
an authoritative or versioned reference.

## mobileOG-db interoperability

mobileOG-db is a curated, actively maintained
database of mobile-genetic-element (MGE) hallmark genes (transposases, integrases, relaxases,
recombinases), searched with DIAMOND. Its hits map onto `functional_features.tsv` with no schema
change: the mobileOG family/gene name to `gene_symbol`, its curated functional role
(integration/excision, replication/recombination/repair, stability/transfer/defense) to
`functional_subcategory`, the transfer/mobility role itself to the existing `mobility_function`
column, `database_name=mobileOG-db`, `database_version` set to the release tag (e.g.
`beatrix-1.6`), and a new documented `functional_category` value: `mobile_genetic_element`.
**Caveat**: mobileOG-db is GPL-3.0 licensed -- that applies to the database file itself if it is
ever bundled with a workflow rather than fetched separately by the user at annotation time.

## Cenote-Taker interoperability

PlasmiCord has no viral/prophage hallmark-gene detection capability today.
Cenote-Taker3 (the successor to the now-deprecated
Cenote-Taker2) is an HMM/homology-based virus and mobile-genetic-element discovery and annotation
pipeline. Its hallmark-gene hits map onto `functional_features.tsv` with no schema change: the
hallmark gene name to `gene_symbol`, its taxonomy call to `product_name`, `database_name`
set to the specific viral-protein database searched, `detection_method=HMM`,
`annotation_engine=Cenote-Taker3`, and a new documented `functional_category` value:
`viral_or_prophage_element` -- kept distinct from mobileOG-db's broader `mobile_genetic_element`
and ISfinder's `insertion_sequence` above, since Cenote-Taker specifically targets viral/phage
hallmark genes, a narrower and different signal from either.

## PlasBench interoperability

Use `plasmicord import plasbench` for native selected-candidate exports, or create
the same manifest directly. The `--mode plasbench` option records the input origin.
See [native adapters](IMPORTS.md) and [shared annotation](../shared_annotation/README.md). PlasBench's
current protein TSV uses zero-based half-open coordinates; convert `start + 1`, retain
`end`, map `sequence_id` to `source_sequence_id`, and map `gene/product/category`
to `gene_symbol/product_name/functional_category`. Do not blindly copy its annotation
digest: the checksum definitions differ. Project-specific recovery scores remain
upstream evidence and are not treated as transmission confidence.
