# Candidate cohorts: Addenbrooke's CRE and UPMC blaNDM-5 outbreak (not yet used for validation)

This is a scoped research note, not an implemented validation cohort. It records what was
confirmed while studying Scherff et al. (the SeqSphere+ real-time plasmid-transmission paper —
see [INPUT_CONTRACT.md's size-correction note](../INPUT_CONTRACT.md#size-corrected-distance-threshold)
for the algorithmic idea PlasmiCord adopted from the same source), so a future contributor does
not have to re-research it from scratch.

## Source

Scherff, Rothgänger, Weniger, Mellmann, Harmsen, "Real-time Plasmid Transmission Detection
Pipeline", *Microbiology Spectrum* (2024), doi:10.1128/spectrum.02100-24 (also posted as bioRxiv
2024.07.09.602722; open access via PMC11619237).

## What's confirmed

Two real, accession-backed cohorts, both used to validate a commercial (Ridom SeqSphere+)
plasmid-transmission module — not independent of that tool's own analysis, but the underlying
sequencing data is public and reusable:

- **Addenbrooke's Hospital carbapenem-resistant cohort**: 81 (of 85 total) carbapenem-resistant
  Gram-negative isolates (*K. pneumoniae*, *E. coli*, *P. putida*, *E. hormaechei*/*roggenkampii*),
  Addenbrooke's Hospital, Cambridge UK, a 6-year collection. **ENA/BioProject PRJEB30134**.
- **UPMC Presbyterian blaNDM-5 outbreak**: 19 isolates, 7 species, Feb 2021–Feb 2023, a US
  hospital outbreak. **BioProject PRJNA981541**. The paper's Discussion uses this cohort's own
  worked example (a transposase insertion making one plasmid 17% larger, raw Mash distance 0.003)
  to justify the size-correction heuristic — the source of PlasmiCord's own
  `--size-correction-per-percent` feature.

## Why these aren't validation cohorts yet

Both cohorts are real and public, but the paper's own comparator (a fixed 0.001 Mash-distance
threshold plus its size correction, and cgMLST-based clonality) is the *method being validated*
in that paper, not an independent ground-truth label PlasmiCord could import directly the way it
imports BLASTn alignment-agreement labels elsewhere (see [CALIBRATION.md](../CALIBRATION.md)).
Using either cohort as a PlasmiCord validation cohort would mean treating the accessions as raw
input and independently re-deriving comparator labels via PlasmiCord's own methodology (e.g.
reciprocal BLASTn agreement, matching the existing 81-plasmid reference cohort's approach) —
not adopting the paper's own alerts as ground truth, since that would be circular for evaluating
a competing/comparable threshold-based method.

Recorded here with confirmed accessions so this doesn't need to be rediscovered later. See
[Klebsiella bloodstream cohort](KLEBSIELLA_BLOODSTREAM_COHORT.md) (the strongest candidate found,
with explicit pair-specific ground truth) and [Shigella cohort](SHIGELLA_COHORT.md).
