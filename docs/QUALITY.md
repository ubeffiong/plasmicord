# Biological evidence and quality rules

The technical gate is followed by a common biological-evidence assessment for every input source. Results appear in `biological_quality.tsv`, the candidate index and the offline report. These research tiers are transparent rules, not calibrated probabilities or a universal plasmid classifier.

Automated checks measure ambiguity, sequence entropy and duplicated terminal sequence. When MOB-typer runs, replicon, relaxase and oriT evidence support plasmid identity. Missing markers leave novel plasmids uncertain rather than rejecting them. A GFA graph can support closure only when a single segment has the exact candidate sequence and a same-orientation self-link. A circular header alone stays tool-reported.

Supply independently produced classification, mapping or contamination evidence with:

```sh
plasmicord run --manifest manifest.tsv --metadata metadata.tsv \
  --quality-evidence quality_evidence.tsv --engine mash --threshold 0.01 --out results_quality
```

The TSV requires `plasmid_id`, `sequence_sha256`, and `evidence_source`. Optional fields are `classification` (plasmid/chromosome/uncertain), `read_breadth` (0–1), `mean_depth` (nonnegative), `chromosome_fraction` (0–1), and `copy_number` (nonnegative). Use an accession, versioned analysis result or reviewed source identifier; preserve the underlying evidence. The candidate digest must match and the evidence file itself is hashed. Missing values remain unmeasured. Merely providing a reads path does not perform mapping.

`copy_number` (e.g. from PCNE's read-coverage-ratio
estimate) is recorded and reported alongside the other evidence, but is **informational only**:
unlike `read_breadth`/`mean_depth`/`chromosome_fraction`, it does not feed the tier policy below.
There is no literature-backed threshold the way there is for those three, and inventing one would
be an unjustified new biological rule -- a very low or very high copy number is not inherently
poor-quality evidence the way low read breadth or high chromosomal contamination is.

PlaScope (a Centrifuge-based contig classifier) is a
plausible upstream source for `classification`: its own output vocabulary is
`plasmid`/`chromosome`/`unclassified`, which maps directly onto this field's vocabulary with
only a rename (`unclassified`→`uncertain`) -- no reinterpretation needed. PlaScope requires a
hand-built, species-specific reference database, so its classification is only as reliable as
that curated reference; treat it like any other supplied evidence source, not as ground truth.

PlasFlow (a composition/k-mer neural-net classifier for
metagenomic contigs) and Plasmer (a random-forest
classifier combining k-mer matching with a replicon-distribution score) are two further
compatible sources: both report a plasmid/chromosome call that maps onto `classification` the
same way. PlasFlow's upstream repository is unmaintained -- treat it as a legacy-compatible
source, not a primary recommendation. Plasmer is actively maintained but requires large
prebuilt reference databases (~32GB RAM), an operational cost worth planning for. Both also
report a numeric per-call probability that `classification` itself cannot carry (it is a plain
label); that probability can optionally be attached via `--external-typing`'s
`predicted_classification_score`/`predicted_classification_call` fields instead -- see
[input contract](INPUT_CONTRACT.md) -- without it influencing this policy either.

PlasForest is a fourth compatible source: a
random-forest classifier using BLAST-homology features (hit count, overlap statistics against a
bundled ~2.5GB reference plasmid set) rather than any of the above tools' k-mer/composition
features, so it can disagree with them in genuinely informative ways. Its categorical output maps
onto `classification` the same zero-reinterpretation way. Like PlasFlow, its upstream repository
has been without a commit for multiple years -- treat it as legacy-compatible, not a primary
recommendation, and budget for its large bundled database.

| Tier | Policy v1.0 |
|---|---|
| High | External plasmid classification, exact graph closure support, breadth ≥0.95, mean depth ≥10, chromosome fraction ≤0.01, and no overriding low-quality flags |
| Moderate | External plasmid classification or MOB plasmid markers, without stronger evidence or overriding warnings |
| Low | Technical low-quality status, low entropy, duplicated terminal sequence, >1% ambiguity, >5% chromosome fraction, breadth <0.95, or mean depth <10 |
| Uncertain | Insufficient positive identity evidence without rejection |
| Rejected | Technical rejection or submitted independent chromosome classification |

Reference accessions establish reference-annotated molecule identity, not read-backed closure. External evidence is validated for identity/schema and retained as supplied; the framework does not independently audit a submitter's classifier or contamination method. High confidence is consequently conditional on that evidence. Long, multi-segment graph traversal, automated read mapping, calibrated species-specific classifiers and empirical validation of these quality tiers remain future work.

`--external-typing` cross-references (MOB-suite cluster IDs, COPLA PTU, host range,
PlasmidFinder Inc-types, sequence-based transmissibility scores, classifier confidence scores,
nearest-PLSDB-reference match -- see [input contract](INPUT_CONTRACT.md)) are reported alongside
these tiers but never feed into this policy: they are opaque external identifiers, not additional
evidence the policy above evaluates.
