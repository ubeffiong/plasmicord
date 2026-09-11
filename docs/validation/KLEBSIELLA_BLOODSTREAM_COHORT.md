# Candidate cohort: Dartmouth-Hitchcock bloodstream K. pneumoniae (not yet used for validation)

This is a scoped research note, not an implemented validation cohort. It records what was
confirmed while studying Long et al./the Nature Communications paper below, so a future
contributor does not have to re-research it from scratch. **This is the strongest validation-cohort
candidate found across all research rounds so far**, because unlike every other candidate
recorded in this directory, it carries an explicit, citable, pair-specific transfer claim rather
than only a population-level rate estimate.

## Source

"Clonal background and routes of plasmid transmission underlie antimicrobial resistance
features of bloodstream *Klebsiella pneumoniae*", *Nature Communications* (2024),
doi:10.1038/s41467-024-51374-x (open access; also on PMC11322185).

## What's confirmed

- **Cohort**: 137 bloodstream *K. pneumoniae* isolates, Dartmouth-Hitchcock Medical Center, New
  Hampshire, USA, January 2017–January 2022. 136 Illumina short-read + 12 ONT long-read hybrid
  assemblies.
- **Accessions**: **BioProject PRJNA1054115** (SRA); isolates also deposited on PathogenWatch;
  individual plasmid assemblies submitted via NCBI BankIT, with accessions listed in the paper's
  Supplementary Data 7 (not independently re-verified row-by-row here — confirm before use).
- **Explicit pair-specific ground truth (the key differentiator)**: the paper's Figure 4
  documents a specific blaCTX-M-15 plasmid transfer history with a named comparator: a
  same-lineage (clonal) ST307→ST307 relationship versus an explicit cross-lineage ST429→ST307
  transfer event, timestamped on a chronological axis (Oct 2019–Nov 2021) with Mash distances
  labeled on the connecting lines in the figure. This is directly analogous to PlasmiCord's own
  cross-cluster-vs-same-cluster discordance distinction — the paper essentially hand-derived the
  same kind of call PlasmiCord automates.
- **Method used by the paper** (for context, not something to reimplement): Mash v2.3 (distance
  ≤0.001) for plasmid clustering, ABRicate/PlasmidFinder for Inc-typing, oriTfinder for
  relaxase/oriT, progressiveMAUVE for structural alignment; RAxML core-genome phylogeny (188,166
  SNPs), BIGSdb MLST/cgMLST, Kleborate for AMR/virulence typing.

## Why it's not a validation cohort yet

PlasmiCord's existing validation cohorts (see [README.md](README.md)) pair real accessions with
an independent, pair-specific comparator label produced by a documented, repeatable method
(BLASTn alignment agreement — see [CALIBRATION.md](../CALIBRATION.md)). This cohort has the
*right shape* of ground truth (an explicit clonal-vs-cross-lineage call), but it exists as a
one-off, hand-curated claim in a single figure/results section, not as a machine-readable label
file with a defined, reproducible method matching PlasmiCord's existing alignment-comparator
methodology.

Using this cohort would require:
1. Obtaining the plasmid/isolate accessions from Supplementary Data 7 and confirming they
   resolve to the specific isolates named in Figure 4.
2. Either treating the paper's own stated clonal/cross-lineage call as the comparator label
   directly (a single labeled pair, or a small handful — check how many such pairs the paper
   actually calls out beyond the headline ST307/ST429 example), or independently re-deriving a
   comparator label via PlasmiCord's own reciprocal-BLASTn methodology across the full 137-isolate
   cohort (a substantially larger validation exercise, closer in scope to the existing 81-plasmid
   reference cohort already documented in [docs/validation/README.md](README.md)).

Either path is a real, separately-scoped follow-up task, not a quick addition — recorded here
with confirmed accessions and the specific pair-level claim so it doesn't need to be
rediscovered later. See [Addenbrooke's/blaNDM-5 cohort](ADDENBROOKES_CRE_COHORT.md) and
[Shigella cohort](SHIGELLA_COHORT.md) for the other candidates found so far.
