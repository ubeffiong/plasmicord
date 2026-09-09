# Population-level summary

`plasmicord population-summary` aggregates outputs a completed `plasmicord run` already
wrote into cohort-level tables. It performs no new inference and requires no new upstream
data -- it is pure aggregation over `plasmid_index.tsv`, `metadata.tsv`,
`plasmid_clusters.tsv`, `biological_quality.tsv` and `functional_features.tsv`.

```sh
plasmicord population-summary --results results_cohort --out-prefix population_summary
```

The target run must have `status: complete` in `run_provenance.json`; each source file is
re-checksummed against the recorded `output_sha256` map before use, so a summary cannot
be silently built from a partially overwritten or stale result directory.

## Outputs

`population_summary.pu_level.tsv` -- one row per plasmid unit:

```text
plasmid_unit n_plasmids n_isolates isolates replicon_types mobility_classes
n_organisms organisms n_locations locations first_date last_date date_range_days
n_resistance_genes resistance_genes n_drug_classes drug_classes interpretation
```

Replicon/mobility values come from `biological_quality.tsv`; resistance genes and drug
classes come from `functional_features.tsv` rows where `headline_eligible == true` (the
same headline-ARG convention already used elsewhere in the report).

`population_summary.metadata_dimension.tsv` -- one row per `(dimension, value)`, for
`dimension` in `organism`, `location`, `chromosomal_cluster` (the three categorical
metadata fields PlasmiCord's own [input contract](INPUT_CONTRACT.md) defines):

```text
dimension value n_isolates n_plasmids n_plasmid_units
n_resistance_genes resistance_genes first_date last_date interpretation
```

## Naming-convention note

This schema is adapted from phac-nml/plasmid_analysis's population-surveillance tables,
which use a dedicated `countries` field. PlasmiCord's metadata contract has no country
column, only a free-text `location`; `location` is reused as-is rather than requiring a
new metadata column. This is an intentional adaptation, not a silently narrowed feature.

Every table entry is a plain aggregation of already-validated run outputs -- it is not a
validated core/accessory-gene claim, an epidemiological linkage claim, or a transmission
inference. See [interpretation](INTERPRETATION.md) for how to read PlasmiCord's outputs
in general.
