# Project review — 8 September 2026

## Assessment

The original project was a small Ubuntu-first proof of concept with useful standalone
clustering, isolate-network and discordance algorithms. It was not yet the independent,
input-flexible functional-cargo framework described in the supplied design notes.
Its final presentation was a short generated Markdown report, TSV files and GraphML
requiring external visualization. There was no interactive HTML report or functional
gene pipeline. The supplied metadata contained placeholders, not a runnable real cohort.

The reviewed version now provides a portable validated core, an offline interactive
report, imported functional-cargo evidence, provenance and regression coverage. It
remains a **research prototype**, not a validated transmission-inference product.

## Correctness findings addressed

| Priority | Original finding | Correction |
|---|---|---|
| High | Unknown chromosome labels could be counted against known labels as cross-cluster events | Exclude unknown labels; report typing coverage |
| High | Missing Mash pairs defaulted to zero and could create false matches | Reject incomplete, invalid and inconsistent matrices |
| High | Contigs were concatenated in extraction/fallback, introducing unsupported joins | Preserve FASTA contigs and block k-mers across boundaries |
| High | Failed/missing assembly extraction silently became missing plasmid evidence | Legacy extraction now stops; fresh output checks prevent cohort mixing |
| Medium | Duplicate IDs and orphan plasmid/isolate mappings could corrupt the graph | Validate mappings, identifier uniqueness and GraphML escaping |
| Medium | Complete-linkage ties depended on input order | Deterministic tie-breaking by member names |
| Medium | Pair–unit counts were described as unique isolate pairs | Report both counts explicitly |
| Medium | Zero crosslinks implied chromosome/plasmid concordance even with missing data | Report no detected crosslinks, with missing-typing limits |
| Medium | Threshold wording implied measured identity and stable unit definitions | Explicit study threshold; sensitivity; run-local PU scope |
| Medium | Single-linkage implied edges without direct distance support | Per-edge minimum distance and direct-threshold flag |

## Added capabilities

- Python CLI and installable package; no third-party Python runtime dependencies.
- Precomputed, assembled-long-read, hybrid and PlasBench-origin manifests through
  one technical gate, without a PlasBench installation or sibling project files.
- FASTA, coordinate, identifier and checksum validation; provenance, duplicate,
  ambiguity, fragmentation and source-quality warnings; rejection records.
- Offline HTML report with interactive network, eligible-ARG filter, gene tracks,
  unit profiles, observed function prevalence, metadata context and downloads.
- Normalized functional features with engine/database provenance and conservative
  headline ARG eligibility. Missing annotation and module completeness stay unresolved.
- Threshold sensitivity, explicit failed/complete run records, input/output checksums,
  deterministic demo and a Windows/Linux CI template.

## Remaining gaps against the supplied design

| Priority | Gap | Present behaviour / next step |
|---|---|---|
| High | Biological plasmid validation | Technical checks only; add classification, contamination, replicon/mobility, read/graph and closure evidence |
| High | Real-cohort validation | Synthetic software tests only; validate reference/plasmid diversity and study-specific thresholds |
| High | Automated annotation | Imports normalized features only; add versioned Bakta/Prokka, a primary AMR caller and replicon/mobility workflow |
| High | Unified internal reconstruction | Legacy MOB shell pipeline remains separate from the new core; route reconstruction through the manifest gate |
| Medium | Specialized native import adapters | Source modes accept the common contract; no automatic PlasBench/MOB/Flye/Unicycler directory discovery |
| Medium | Raw-read preparation | Not implemented; reads require upstream assembly, polishing and plasmid extraction |
| Medium | Shared annotation/viewer package with PlasBench | Contract documented; no shared versioned engine/viewer package extracted yet |
| Medium | Annotation evaluation denominators | No per-plasmid caller-completion matrix; observations cannot establish gene absence or validated core/accessory calls |
| Medium | Rich functional/network controls | Initial metadata/ARG/category/search filters present; dedicated drug-class, mechanism, mobility, date-range and module controls remain |
| Medium | Orthology and metabolic completeness | Supplied terms preserved; eggNOG/KO/domain calls and defined-component completeness not computed |
| Medium | Annotation caching and reconciliation | No sequence/tool/database/parameter cache or cross-database allele reconciliation |
| Medium | Pair confirmation and epidemiological modelling | No alignments, synteny, containment, direction, uncertainty model or risk score |
| Medium | Scale and orchestration | Dense matrices and pair enumeration; no cohort-scale performance validation or workflow manager |
| Low | Stable cross-run unit registry | Deterministic IDs for a fixed cohort only; registry/version matching remains future work |

PlasBench was reviewed but not changed. Its reconstruction-recovery view remains a
separate responsibility. Project 2 is a continuation in scientific scope, not an
application dependency. See [POSITIONING.md](POSITIONING.md).

## Validation record

- Original clustering and network/discordance test suites passed before and after changes.
- Ten new regression tests passed: unknown clusters, missing/invalid distances,
  deterministic complete-linkage ties, FASTA/coordinate/checksum validation, loose
  ARG exclusion, contig-boundary handling, empty/singleton runs and HTML escaping.
- Windows Python 3.14: full exact-k-mer synthetic run passed.
- Ubuntu/WSL system Python: full Mash synthetic run passed, invoking the existing
  Mash executable from an available bioinformatics environment. No PlasBench Python
  package or workflow was invoked. This checks Mash integration, not reconstruction.
- Both engines produced 5 isolates, 5 plasmids, 2 units, 4 sharing pairs and 3 unique
  cross-cluster pairs. Annotations were illustrative synthetic fixtures.
- Python wheel built, installed to an isolated target and ran its full demo with
  the offline report template included.
- The HTML opened offline in Chromium; eligible-ARG and cross-cluster filters,
  unit selection, gene details and edge details were exercised without browser errors.

No real-data MOB reconstruction, biological annotation engine or surveillance cohort
was run in this review. The GitHub login could create the public repository but lacked
the `workflow` scope required to upload an active Actions workflow. The tested project
is published with [ci-template.yml](ci-template.yml); CI is not enabled. A maintainer
with workflow write permission can copy it to `.github/workflows/ci.yml`.
