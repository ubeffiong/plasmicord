# Public-reference validation — 8 September 2026

These results use real, accession-backed sequences. The separate seeded demo remains synthetic. PlasmiCord copied references into an independent input package and did not invoke or import the PlasBench application.

## Cohort and main run

[Reference metadata](reference_metadata.tsv) lists 32 isolates, their organisms, BioProjects and assembly accessions. [Reference manifest](reference_manifest.tsv) lists 81 NCBI-designated plasmids and sequence digests. All 32 locally available cohort samples were retained; [source records](reference_sources.json) preserve reference FASTA and NCBI sequence-report hashes. FASTAs are not bundled here; manifest paths describe the original standalone input package. Obtain corresponding sequences by their versioned NCBI accessions.

The core run used Mash 2.3, k=21, sketch size 10,000, complete linkage and distance 0.01. It produced 63 units, including 49 singletons, and nine isolate-sharing edges. NCBI molecule annotations supplied plasmid-classification evidence. No independent read mapping, contamination measurements or graph closure were supplied.

All 32 isolates lacked comparable chromosome-cluster labels. Cross-cluster inference was therefore not evaluable; zero reported crosslinks do not imply concordance. Dates and locations were not invented.

## Reference-relatedness calibration

The BioProject partition was assigned deterministically before scoring. Cross-split exact sequence duplicates were excluded. Only different-isolate pairs within a BioProject were compared. BLASTn 2.17.0 reciprocal alignments provided comparator labels, using the rules in [CALIBRATION.md](../CALIBRATION.md). The labels are independent of Mash calculations but are not epidemiological transmission truth.

| Partition | BioProjects | Positive | Negative | Uncertain |
|---|---:|---:|---:|---:|
| Training | 4 | 5 | 149 | See pair-level labels |
| Held out | 3 | 9 | 111 | See pair-level labels |

There were 289 pairs: 274 definite and 15 uncertain. Training selected **0.001** from the prespecified sweep 0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05. The holdout had TP=9, TN=111, FP=0, FN=0. The very small positive sample and within-study dependence prevent claims of universal accuracy. Several training cutoffs tied; the selection rule chose the smallest. Do not deploy this cutoff universally or treat it as a transmission probability.

The [calibration report](calibration/CALIBRATION.html), [metrics](calibration/calibration_metrics.tsv), [pair labels](labels.tsv), [label provenance](label_provenance.json), and compressed [raw BLAST output](blastn.tsv.gz) are included. The gzip contains the exact bytes identified by the raw-output SHA-256.

Replay threshold fitting from this repository, without bioinformatics tools or sequence downloads:

```sh
python plasmicord.py calibrate --matrix docs/validation/plasmid_matrix.tsv \
  --index docs/validation/plasmid_index.tsv --labels docs/validation/labels.tsv \
  --cohort-kind real --engine mash --out results_calibration_replay
```

The portable index retains IDs and sequence hashes; machine-specific paths were removed. The [source run record](run_provenance.json) preserves original run/index hashes and checksums of these portable calibration inputs.

## Real automatic annotation exercise

The essential profile ran on [NZ_CP054237.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP054237.1) (131,660 bp) and [NZ_CP054238.1](https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP054238.1) (1,546 bp), from assembly GCF_013372425.1 / BioProject PRJNA636382.

- Prokka 1.15.6 with its bundled database produced 151 CDS on the larger plasmid and zero on the smaller.
- AMRFinderPlus 4.2.7, database 2026-08-07.1, produced 12 functional matches on the larger plasmid, including seven headline-eligible ARGs: blaTEM-1, sul2, aph(3'')-Ib, aph(6)-Id, dfrA1, aadA1 and sul1. The remaining matches concern biocide/metal stress. The smaller plasmid had zero AMRFinder matches.
- MOB-typer 3.1.9 reported IncFIA/IncFIB/IncFII/IncQ1 replicons and a conjugative prediction on the larger plasmid, and Col(MG828) with a mobilizable prediction on the smaller. These are tool predictions, not demonstrated transfer phenotypes.
- All six stages completed. A repeat run reused all six cached stages. Both candidates received moderate marker-supported quality; neither was called read-validated or confirmed closed.

[Normalized features](annotation/functional_features.tsv), [completion records](annotation/annotation_status.tsv), [provenance](annotation/provenance.json), [typing](annotation/mobility_typing.json), [biological evidence](annotation/biological_quality.tsv), and [protein export](annotation/plasbench_proteins.tsv) are included. Published provenance omits local executable/database/output paths while retaining tool versions and content hashes. The MOB database's release label was unavailable: its initialization date was 3 September 2026 and its full content digest is recorded. No precise release number is invented.

This is an integration test on two real plasmids, not an annotation-accuracy benchmark across all 81 references. Zero detected CDS/ARGs means that the specified tool completed with zero calls, not that biology guarantees absence. Bakta's database download failed checksum/transport checks; Bakta is supported by the implementation and parser tests but was not run on these references.

## Software checks

Ten standalone regressions, 17 gap-workflow tests, both original algorithm suites, Windows synthetic runs, WSL Mash/annotation runs, native PlasBench-export import, independent package builds/install checks and offline Chromium report checks were performed. PlasBench native import used its existing synthetic audit export and is not counted as real-cohort evidence.
