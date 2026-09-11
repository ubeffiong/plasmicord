# Candidate cohort: MDU PHL Shigella (not yet used for validation)

This is a scoped research note, not an implemented validation cohort. It records what was
confirmed while studying the supplementary code for Müller, Wick, Judd, Williamson, Bedford,
Howden, Duchêne, Ingle, *Quantifying plasmid movement in drug-resistant Shigella species using
phylodynamic inference*, PLOS Pathogens 21(12): e1013621, 2025 (PMC12677775), so a future
contributor does not have to re-research it from scratch.

## What's confirmed

- **Cohort**: 789 *Shigella sonnei* and 297 *S. flexneri* isolates (excluding serotype 6),
  January 2016–December 2020, received by the Microbiological Diagnostic Unit Public Health
  Laboratory (MDU PHL), Doherty Institute, Victoria, Australia.
- **Accessions — public today, not on request**: the paper's Data Availability Statement
  points to supplementary **S1 Table** (sonnei) and **S2 Table** (flexneri), each a CSV with
  `run_accession`, `sample_accession`, `study_accession` per isolate. All sampled rows share
  **BioProject PRJNA857526** ("Whole genome sequencing of Escherichia/Shigella sp as part of
  MDU-PHL surveillance activities in Victoria, Australia" — 1,743 BioSamples, 4,101 SRA runs,
  verified live on NCBI). S2 Table (flexneri) accessions were not independently re-verified
  row-by-row — confirm before use.
- **Repo contents**: `Applications/Shigella/` in the paper's supplementary code has plasmid-presence/date
  matrices (e.g. `sonnei_plasmid_info.tsv`: isolate ID, date, binary presence per reference
  plasmid) — useful structure, but **no accession numbers**; those live only in the paper's
  S1/S2 tables. `xmls/`/`supplementalxmls/` hold ~90+ BEAST2 XML configs encoding taxon-date
  pairs.

## Why it's not a validation cohort yet

PlasmiCord's existing validation cohorts (see [README.md](README.md)) pair real accessions
with an independent, pair-specific comparator label (BLASTn alignment agreement). This
paper's output is **population-level transfer-rate estimates** (events/year, per plasmid),
not a citable table of discrete transfer events with source/destination isolate pairs and
dates — the "events" live inside BEAST2 posterior tree-space, not as an exported label file.

Using this cohort for PlasmiCord validation would need one of:
1. Treating the repo's own plasmid-presence/date TSVs as a **weaker same-lineage-presence
   proxy** label (not a transfer-event label) — the more tractable option.
2. Rerunning the published BEAST2/CoalPT XMLs to extract posterior transfer events directly —
   a substantial undertaking (BEAST2/MCMC, per-lineage tree inference), out of PlasmiCord's
   stdlib-only scope for the tool itself, though the *output* of such a rerun could in
   principle be imported as comparator labels the same way existing alignment labels are.

Either path is a real, separately-scoped follow-up task, not a quick addition — recorded
here with confirmed accessions so it doesn't need to be rediscovered later. See
[CALIBRATION.md](../CALIBRATION.md#relationship-to-phylodynamic-transfer-rate-inference) for
how this relates to PlasmiCord's existing threshold-calibration approach.
