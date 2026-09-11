# Study-specific threshold calibration

Calibration selects a distance cutoff against independently supplied labels. It never uses existing Mash clusters as truth.

```sh
plasmicord calibrate --matrix results_cohort/plasmid_matrix.tsv \
  --index results_cohort/plasmid_index.tsv --labels labels.tsv \
  --cohort-kind real --engine mash --linkage complete \
  --thresholds 0,0.001,0.002,0.005,0.01,0.02,0.05 --out calibration
```

Label TSV columns: `plasmid_a,plasmid_b,label,split,group_id,evidence_source,evidence_type`.
Labels are positive/negative/uncertain; splits are train/validation. Evidence types are reviewed_plasmid_sharing, experimental, independent_alignment or synthetic_truth. Synthetic truth cannot be presented as a real cohort.

Allocate entire independent studies/patients to one split before fitting. Candidate IDs, isolates, exact sequence digests and group IDs cannot overlap splits. The default requires both classes and at least three groups in each split. Insufficient evidence yields **no selected threshold**. Reducing that guardrail is possible for exploratory work but does not establish adequate validation.

For each training threshold, candidates are clustered within the training partition and pair co-membership is scored. Selection maximizes training balanced accuracy, then F1, then prefers the smaller threshold. Only the selected threshold is evaluated on the held-out partition. Uncertain pairs are excluded from scoring; their candidates remain in the partition's clustering. No validation metric participates in threshold selection.

A neighbouring `run_provenance.json` verifies source completion, engine and matrix checksum and carries sketch parameters into the calibration record. Without one, the engine is marked submitter-declared. Calibration only applies to the recorded engine, sketch parameters, linkage and population. It is not a direct-transmission probability.

The final `CALIBRATION.html` presents the selection, limitations, training sweep and holdout confusion counts. TSV metrics, labels and JSON provenance accompany it. See [public validation](validation/README.md).

The independent reference workflow is:
1. `scripts/prepare_reference_validation.py`: extract NCBI-designated plasmids from a local accession-backed assembly package.
2. Run PlasmiCord to calculate distances.
3. `scripts/label_reference_alignments.py`: whole-BioProject deterministic split, exclusion of cross-split exact duplicates, reciprocal BLASTn alignment comparator labels.
4. Run the calibration command.

Alignment-positive defaults require ≥99.5% weighted HSP identity and ≥90% union coverage of both sequences in both directions. Negative requires <50% coverage in both directions; intermediate cases are uncertain. Repetitive HSPs are not a one-to-one synteny assessment. These comparators evaluate sequence-relatedness agreement; they are not reviewed epidemiological transmission labels. Pair observations within studies are correlated; no population confidence interval is claimed.

## Relationship to phylodynamic transfer-rate inference

This calibration fits a single static distance threshold against labeled pairs -- a
lightweight, cohort-agnostic heuristic. It is not a substitute for full phylodynamic
inference. Müller et al.'s joint coalescent + plasmid-transfer model (BEAST2/CoalPT; see
their *Quantifying plasmid
movement in drug-resistant Shigella species using phylodynamic inference*, PLOS Pathogens
2025) jointly models chromosome and plasmid trees to estimate transfer/loss rates over time
-- structurally more rigorous, but requiring BEAST2/MCMC and per-lineage tree inference, out
of PlasmiCord's stdlib-only scope. See [the Shigella cohort note](validation/SHIGELLA_COHORT.md)
for a real, accession-backed dataset from the same paper recorded as a scoped future
validation-cohort candidate.
