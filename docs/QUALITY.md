# Biological evidence and quality rules

The technical gate is followed by a common biological-evidence assessment for every input source. Results appear in `biological_quality.tsv`, the candidate index and the offline report. These research tiers are transparent rules, not calibrated probabilities or a universal plasmid classifier.

Automated checks measure ambiguity, sequence entropy and duplicated terminal sequence. When MOB-typer runs, replicon, relaxase and oriT evidence support plasmid identity. Missing markers leave novel plasmids uncertain rather than rejecting them. A GFA graph can support closure only when a single segment has the exact candidate sequence and a same-orientation self-link. A circular header alone stays tool-reported.

Supply independently produced classification, mapping or contamination evidence with:

```sh
plasmicord run --manifest manifest.tsv --metadata metadata.tsv \
  --quality-evidence quality_evidence.tsv --engine mash --threshold 0.01 --out results_quality
```

The TSV requires `plasmid_id`, `sequence_sha256`, and `evidence_source`. Optional fields are `classification` (plasmid/chromosome/uncertain), `read_breadth` (0–1), `mean_depth` (nonnegative), and `chromosome_fraction` (0–1). Use an accession, versioned analysis result or reviewed source identifier; preserve the underlying evidence. The candidate digest must match and the evidence file itself is hashed. Missing values remain unmeasured. Merely providing a reads path does not perform mapping.

| Tier | Policy v1.0 |
|---|---|
| High | External plasmid classification, exact graph closure support, breadth ≥0.95, mean depth ≥10, chromosome fraction ≤0.01, and no overriding low-quality flags |
| Moderate | External plasmid classification or MOB plasmid markers, without stronger evidence or overriding warnings |
| Low | Technical low-quality status, low entropy, duplicated terminal sequence, >1% ambiguity, >5% chromosome fraction, breadth <0.95, or mean depth <10 |
| Uncertain | Insufficient positive identity evidence without rejection |
| Rejected | Technical rejection or submitted independent chromosome classification |

Reference accessions establish reference-annotated molecule identity, not read-backed closure. External evidence is validated for identity/schema and retained as supplied; the framework does not independently audit a submitter's classifier or contamination method. High confidence is consequently conditional on that evidence. Long, multi-segment graph traversal, automated read mapping, calibrated species-specific classifiers and empirical validation of these quality tiers remain future work.

`--external-typing` cross-references (MOB-suite cluster IDs, COPLA PTU, host range -- see
[input contract](INPUT_CONTRACT.md)) are reported alongside these tiers but never feed
into this policy: they are opaque external identifiers, not additional evidence the
policy above evaluates.
