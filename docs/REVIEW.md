# PlasmiCord review — v0.4.0, 8 September 2026

PlasmiCord now includes implementations addressing the five requested gaps. It remains a research framework; neither a sharing edge nor a quality tier is a validated transmission inference.

## Requested gap fixes

| Gap | Implementation | Scope |
|---|---|---|
| Automatic annotation | Essential Prokka/Bakta + AMRFinderPlus + MOB-typer profile; explicit subsets; normalized CDS/AMR evidence, completion records and content-keyed cache | Real Prokka/AMRFinder/MOB integration exercise. Bakta execution remains unverified locally after its database download failed |
| Biological quality | Complexity/ambiguity/terminal-repeat checks; MOB markers; exact GFA closure; checksum-linked external classification, read and contamination evidence | NCBI references exercise identity evidence. Quality tiers remain heuristic research rules |
| Native imports | MOB-recon, Flye, Unicycler, selected PlasBench and generic FASTA | Common gate, grouping/selection safeguards, source hashes and circularity provenance |
| Shared annotation | Independently buildable plasmid-annotation-core wheel/CLI; shared engine used by PlasmiCord; PlasBench protein export | No application dependency. PlasBench has not adopted the library internally; its viewer stays separate |
| Real-cohort calibration | Training-only selection, group/isolate/exact-sequence leakage guards, held-out scoring and abstention | 81 accession-backed plasmids from 32 isolates; independent BLAST comparisons across seven BioProjects. Sequence-relatedness agreement, not transmission truth |

## Detailed final report added in v0.4.0

The new dashboard adds deterministic, evidence-linked module findings; dynamic selected-isolate,
network and candidate interpretations; colour-coded statuses; chart exports; a distance heatmap;
complete per-isolate reports; sortable/searchable/paginated tables; optional matching calibration;
and a recursive output tree with safe previews, metadata, checksums and original downloads.
A portable ZIP includes the result directory. Diagnostic reports distinguish unfinished runs from
actual rejection or completed zero-hit searches. See [REPORT_GUIDE.md](REPORT_GUIDE.md).

The presentation was generated and checked on the existing synthetic demo, 81-plasmid reference
cohort and two real annotated plasmids. Analysis was not repeated merely to refresh HTML.
Seven new report regressions cover integrity, failure states, coverage, calibration and packaging;
browser assertions exercise the visible controls, evidence links, charts and print pagination.

## Original correctness fixes retained

Unknown chromosome labels are excluded from cross-cluster calls. Incomplete or inconsistent Mash matrices fail instead of creating zero-distance matches. Contig boundaries are preserved; identifier, coordinate and checksum errors fail early. Complete-linkage ties are deterministic. The report separates unique isolate pairs from pair–unit observations and shows direct distance support. Missing annotations and typing remain explicit. Fresh output directories and failure provenance prevent stale-result mixing.

## Outputs and presentation

Open **REPORT.html** from a completed run. It is a self-contained offline report with interactive isolate network, metadata/ARG filters, edge evidence, chromosome-cluster discordance, candidate gene tracks, function prevalence, biological-quality evidence and annotation completion. Supporting TSV, JSON and GraphML downloads retain the evidence. REPORT.md is the plain-text counterpart.

Automatic runs retain native tool outputs/logs, annotation provenance and a PlasBench-compatible protein export. Separate **CALIBRATION.html** presents the training sweep, chosen threshold, held-out confusion counts and provenance; it never silently changes an analysis threshold. See [README](../README.md), [annotation](ANNOTATION.md), [quality](QUALITY.md), [imports](IMPORTS.md), and [calibration](CALIBRATION.md).

## Validation

The [public validation record](validation/README.md) contains accessions, comparator labels, parameters and limitations. The reference analysis at Mash distance 0.01 with complete linkage produced 63 units and nine isolate-sharing edges. All 32 isolates lacked supplied chromosome clusters, so chromosome discordance was not evaluable. No chromosome labels, dates or locations were fabricated.

Calibration selected 0.001 using four training BioProjects, then scored 9 true-positive and 111 true-negative pairs in three held-out BioProjects, with zero errors against the specified alignment comparator. Only five training positives and nine holdout positives were available; 15 intermediate pairs were excluded from scoring. Perfect agreement on this small exercise does not establish universal accuracy.

Software validation covers ten standalone regressions, the new gap-workflow tests, both original algorithm suites, synthetic end-to-end reports, package builds and isolated installs. Synthetic annotations remain illustrative fixtures; real annotation evidence is separately identified in the validation record.

## Remaining work

- Independently reviewed transmission labels, comparable chromosome typing and external-population validation.
- Automated raw-read assembly/polishing, mapping and contamination classification. Current quality assessment consumes external measurements and automatically evaluates sequence, marker and simple graph evidence.
- Complex graph closure, synteny and structural confirmation remain future work. A new length/similarity heuristic (`containment_candidates.tsv`) flags candidate nested-plasmid pairs, but Mash/k-mer distance is a poor proxy for size-disparate containment (a small plasmid nearly wholly nested in a much larger one), so it only catches near-equal-length, high-similarity pairs; alignment-based (MUMmer4-style) confirmation is unimplemented and out of stdlib-only scope. A replicon does not prove a complete plasmid.
- Empirically calibrated quality tiers, phenotypic AMR validation and organism-specific mutation analysis.
- Orthology, validated conserved/accessory functions and defined-component module completeness.
- Shared viewer extraction and explicit adoption of the independent library inside PlasBench.
- Dense-matrix scaling, cohort-scale performance testing, a stable cross-run unit registry and orchestration.
- The older MOB shell pipeline remains separate; its outputs can now enter the common gate through the native importer.
- External typing/taxonomy cross-references (MOB-suite cluster IDs, COPLA PTU, host-range fields) can be attached via a checksum-linked `--external-typing` TSV; these remain opaque external identifiers that never influence PlasmiCord's own quality/confidence tiers.
- Threshold-margin diagnostics (`network.edge_evidence.tsv`, `plasmid_clusters.tsv`) expose how close a link or plasmid unit is to its distance threshold, surfacing single-linkage chaining risk without changing the clustering algorithm itself.
- An optional `plasmicord population-summary` step aggregates existing run outputs into cohort-level tables (plasmid-unit and metadata-dimension summaries), adapted to PlasmiCord's own metadata contract (`location` in place of a dedicated country field). See [PlasBench recommendations](PLASBENCH_RECOMMENDATIONS.md) for lessons from the same competitor research that apply to PlasBench's reconstruction-benchmarking mission instead.

GitHub Actions remains disabled because the available authorization lacks workflow-write scope. The [CI template](ci-template.yml) includes the new tests and can be installed by a maintainer with that permission.
