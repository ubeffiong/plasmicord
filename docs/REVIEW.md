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
- Complex graph closure, synteny and structural confirmation remain future work. A replicon does not prove a complete plasmid.
- A new length/similarity heuristic (`containment_candidates.tsv`, always computed) flags candidate nested-plasmid pairs, but Mash/k-mer distance is a poor proxy for size-disparate containment (a small plasmid nearly wholly nested in a much larger one), so it only catches near-equal-length, high-similarity pairs; alignment-based (MUMmer4-style) confirmation is unimplemented and out of stdlib-only scope. Widening `--containment-max-ratio` to 1.0 would also flag equal-length pairs; the default 0.95 excludes them deliberately.
- Empirically calibrated quality tiers, phenotypic AMR validation and organism-specific mutation analysis.
- Orthology, validated conserved/accessory functions and defined-component module completeness.
- Shared viewer extraction and explicit adoption of the independent library inside PlasBench.
- Dense-matrix scaling, cohort-scale performance testing, a stable cross-run unit registry and orchestration.
- The older MOB shell pipeline remains separate; its outputs can now enter the common gate through the native importer.
- External typing/taxonomy cross-references (MOB-suite cluster IDs, COPLA PTU, host-range fields) can be attached via a checksum-linked `--external-typing` TSV; these remain opaque external identifiers that never influence PlasmiCord's own quality/confidence tiers.
- Threshold-margin diagnostics (`network.edge_evidence.tsv`, `plasmid_clusters.tsv`) expose how close a link or plasmid unit is to its distance threshold, surfacing single-linkage chaining risk without changing the clustering algorithm itself.
- An optional `plasmicord population-summary` step aggregates existing run outputs into cohort-level tables (plasmid-unit and metadata-dimension summaries), adapted to PlasmiCord's own metadata contract (`location` in place of a dedicated country field). See [PlasBench recommendations](PLASBENCH_RECOMMENDATIONS.md) for lessons from the same competitor research that apply to PlasBench's reconstruction-benchmarking mission instead.
- A second round of repo study (PlasTrans, PlasAnn, Ecological-Complexity-Lab's Plasmid_multilayer_networks, TaDReP, plascad, Plasmids-material) found none duplicate PlasmiCord's core mission. Adopted: `--external-typing` now also accepts a generic sequence-based transmissibility score/call (e.g. from PlasTrans's codon-usage CNN); a new `tadrep` native importer adapter (schema verified against TaDReP's own source, not just its README -- its own `plasmids.info.tsv` has a confirmed upstream key-mismatch bug in its `CDS` column, not imported here but worth knowing if a future contributor wires that file in too); and an opt-in `--multilayer-network` export (one edge per plasmid unit, tagged `same_cluster`/`cross_cluster`/`unknown_cluster`), applying Ecological-Complexity-Lab's "layer = categorical partition" pattern to PlasmiCord's own isolate/cluster data rather than their plasmid-similarity/host-layer model (which doesn't map onto PlasmiCord's node type) -- this is edge-relabeling by existing fields, not a true node-multiplex/projection framework like `pymnet`/`multinet`. PlasAnn's oriT/oriV/transposon calls map onto the existing functional-feature schema with no code change (see [input contract](INPUT_CONTRACT.md)).
- Not adopted: plascad's MOB/MPF-HMM classification + SARG ARG detection is redundant with the already-imported MOB-typer output and already-wrapped AMRFinderPlus. Ecological-Complexity-Lab's ~2,000-plasmid dairy-cow plasmidome dataset (real, MIT-licensed) is domain-mismatched for PlasmiCord's clinical-isolate/outbreak-surveillance validation use case (host-microbiome ecology, not AMR surveillance).
- Müller et al.'s BEAST2/CoalPT joint coalescent + plasmid-transfer inference (see Plasmids-material, PLOS Pathogens 2025) is a structurally more rigorous, complementary alternative to `calibration.py`'s static distance-threshold fit; reimplementing BEAST2/MCMC is out of stdlib-only scope. The same paper's real, public, accession-backed Shigella cohort (BioProject PRJNA857526) is recorded as a scoped future validation-cohort candidate — see [docs/validation/SHIGELLA_COHORT.md](validation/SHIGELLA_COHORT.md).
- A third round of study (Scherff et al./SeqSphere+, wtmatlock/plasmid-network-analysis, PlaScope, PlasmidFinder, harmfull_plasmids, pentamorfico/plasmid_network, a Nature Communications *K. pneumoniae* paper, an MDPI mcr-gene paper) again found no duplication of PlasmiCord's core mission. Adopted: an opt-in `--size-correction-per-percent` clustering-threshold loosener for plasmid pairs with large length differences, formula verified against Scherff et al.'s own worked example after an automated first read produced a self-contradictory summary (see [input contract](INPUT_CONTRACT.md#size-corrected-distance-threshold)); `--external-typing` now also accepts PlasmidFinder Inc-type calls (independently-sourced from MOB-typer's own database, so disagreement is evidence, not something reconciled); a documented mapping from PlaScope's plasmid/chromosome/unclassified output onto `--quality-evidence`'s existing classification field (see [QUALITY.md](QUALITY.md)); and four dashboard visualizations inspired by the alignment-track, network, and encoding conventions in this round's sources — a pairwise gene-track comparison for a selected network edge, an isolate timeline with distance-labeled sharing edges, a log/linear toggle for the unit-size chart, and node-size dual-encoding for the network view (all vanilla-JS/inline-SVG, no new dependency; see `python/report.js`/`report_dashboard.js`).
- Not adopted: harmfull_plasmids is agent-based population-dynamics *simulation* code (no genomic input mode at all) — a name collision with "harmful plasmid" sequence classifiers, not an applicable `--external-typing` source despite the superficially matching name. wtmatlock/plasmid-network-analysis's environmental/livestock F-type plasmid dataset is domain-mismatched for clinical outbreak validation (consistent with the earlier dairy-cow-dataset precedent) and carries no license. The MDPI mcr-gene paper's 5,549-plasmid, 85-country dataset is a real but aggregate public-database meta-analysis, not a single controlled cohort with alignment-based comparator labels — lower priority than the two cohorts below; its bubble-matrix (category×category, size-encoded) visual is noted as a lower-priority future population-summary chart idea, not built this round.
- Two further real, accession-backed cohorts are recorded as scoped future validation-cohort candidates: [docs/validation/KLEBSIELLA_BLOODSTREAM_COHORT.md](validation/KLEBSIELLA_BLOODSTREAM_COHORT.md) (BioProject PRJNA1054115 — the strongest candidate found across all rounds, with an explicit pair-specific clonal-vs-cross-lineage transfer claim, not just a population-level rate) and [docs/validation/ADDENBROOKES_CRE_COHORT.md](validation/ADDENBROOKES_CRE_COHORT.md) (BioProjects PRJEB30134 and PRJNA981541, from the same paper as the size-correction feature).

- A fourth round of study (seqviz, plasmidID, PlasFlow, ISfinder-sequences, Recycler, mobileOG-db, Plasmer, motif, SpliceCraft, tiptoft) again found no duplication of PlasmiCord's core mission. Adopted: a circular gene-map view in the dashboard (backbone tick ring plus one annotation ring, toggled alongside the existing linear track) reusing seqviz's SVG arc-math *technique* only -- the library itself requires React and can't be used in a zero-dependency single-file report -- combined with plasmidID's ring-layout template (its coverage/contig-alignment rings were skipped: PlasmiCord has no per-base coverage data to draw them); `--external-typing` now also accepts a generic `predicted_classification_score`/`predicted_classification_call` pair for composition-based plasmid/chromosome classifiers such as PlasFlow and Plasmer (both otherwise map their categorical call onto `--quality-evidence`'s classification field the same way PlaScope already does -- see [QUALITY.md](QUALITY.md) -- this new field is a distinct evidentiary axis from `predicted_transmissibility_score`, not a reuse of it); and two further docs-only `functional_features.tsv` mapping guides for ISfinder-sequences (`functional_category=insertion_sequence`; the mirror itself carries no license and has been stale since 2020, documented as a convenience source only) and mobileOG-db (`functional_category=mobile_genetic_element`; GPL-3.0 applies if the database file is ever bundled) -- see [input contract](INPUT_CONTRACT.md).
- Not adopted: motif and SpliceCraft are both wet-lab/design-stage sequence workbenches (cloning, primer design, MSA editing) despite ambiguous names suggesting otherwise -- `motif` is not a motif-finder -- and are not applicable to a surveillance/reporting framework. Recycler is a pure assembly-graph-internals reconstruction algorithm (SPAdes FASTG + BAM coverage-consistency cycle detection); it confirms rather than changes PlasmiCord's existing input-boundary framing (candidates are assumed already reconstructed). tiptoft predicts plasmid replicons directly from raw long reads before assembly -- out of scope for a tool that starts only after reconstruction, though worth knowing as a diagnostic hint: a plasmid an assembler dropped might still show up in a tiptoft scan of the same isolate's raw reads.

- A fifth round of study (plsdb, repp [Lattice-Automation and jjti confirmed the same project, development moved from the original author's repo to his company's], PlasForest, roundabout, Cenote-Taker2→3, octopus3, yapv, PCNE, Plasmid-Planner) again found no duplication of PlasmiCord's core mission. Adopted: `--external-typing` now also accepts a nearest-reference cross-reference against PLSDB (~72,360 curated plasmids, CC-BY) -- `plsdb_nearest_accession`/`plsdb_nearest_distance`/`plsdb_nearest_host`, populated by running `mash dist` locally against PLSDB's own downloadable Mash sketch, no new dependency; PlasForest is documented as a fourth compatible `--quality-evidence` classifier alongside PlaScope/PlasFlow/Plasmer (a genuinely different BLAST-homology feature set, though also legacy-compatible/unmaintained -- see [QUALITY.md](QUALITY.md)); `--quality-evidence` gains an optional `copy_number` field (e.g. from PCNE) that is recorded and reported but deliberately never influences quality tiers, since no literature-backed threshold exists for it the way one does for read breadth/depth/contamination; and a docs-only `functional_features.tsv` mapping guide for Cenote-Taker3's viral/prophage hallmark-gene hits (`functional_category=viral_or_prophage_element`, kept distinct from mobileOG-db's broader MGE category and ISfinder's insertion-sequence category). While implementing this round's `copy_number` validation, a pre-existing gap was also fixed: `load_evidence()`'s `mean_depth` check (and the new `copy_number` check that mirrored it) crashed with a raw, unhelpful Python `ValueError` on non-numeric input instead of the clean message every other numeric field in this codebase already gives; both now fail cleanly.
- Not adopted: repp and Plasmid-Planner are wet-lab DNA-assembly/cloning *construction-planning* tools (desired sequence in, build protocol out), the same category as SpliceCraft/motif already rejected. roundabout is a real, actively-maintained plasmid-outbreak-clustering pipeline, but its annotation stack (AMRFinderPlus, PlasmidFinder) is already directly integrated into PlasmiCord -- this round confirms rather than extends that integration; its plotting stack (PyGenomeViz/MinkeMap/DaisyBlast) renders separate images, not a portable inline-SVG technique. octopus3's reference-concordance colony-QC idea (map reads to expected sequence, flag variants) is conceptually adjacent but too generic to turn into a concrete new evidence field without inventing an under-specified rule -- noted only as a framing idea for possible future work. yapv is a genuinely dependency-free (no React/build-step, unlike seqviz) plasmid-map renderer -- technically inlinable without violating PlasmiCord's single-file constraint -- but 4.5 years without a commit, and offers no visual capability (GC-content ring, restriction sites) beyond PlasmiCord's own circular-map view already shipped in the prior round; its generic Track/Marker/leader-line schema is noted as a possible future internal convention if a second annotation-ring type is ever added, not built now.

## v0.5.0, 11 September 2026: production-readiness audit, comprehensive test suite and dashboard completion

A full production-readiness audit found three "silo" outputs -- `containment_candidates.tsv`,
`typing_crossreference.tsv` and `network.multilayer_edges.tsv`/`.graphml` -- were computed and
written to disk but never threaded into `report_data.json`, making them invisible to anyone
using the dashboard rather than reading raw TSVs directly. Fixed: `report.py::build_report()`
now accepts and embeds all three (with a defensive disk-read fallback for the `plasmicord
report` regeneration path), `report_evidence.py` gained matching `describe()` catalog entries
and `finding()` rules, and the dashboard's Sharing and Quality sections now have dedicated
sortable/searchable/exportable tables for containment candidates, multilayer network edges and
the external typing cross-reference, not just a findings-grid summary sentence. Also fixed in
the same audit: an unguarded `float()` in `calibration.py`'s `--thresholds` parsing, a stale
gap-workflow test count in the validation record, and a missing visual cue disabling the
now-inapplicable "Track zoom" control in circular gene-map view.

Test coverage was substantially widened: `test_cli_end_to_end.py` invokes `plasmicord.py` as a
real subprocess for every subcommand (argparse parsing, exit codes and stderr text included --
previously every test called internal functions directly with a hand-built `argparse.Namespace`,
never exercising real CLI behavior), and `test_edge_cases.py` covers malformed manifests,
dangling file references, duplicate IDs, non-UTF-8 input, an all-rejected cohort, a missing
`mash` binary, out-of-order CLI bounds and a tampered `run_provenance.json`. A new
`scripts/full_feature_demo.py` (`make full-demo`) exercises every dashboard capability in one
synthetic run. A new opt-in, network-touching `scripts/fetch_real_cohort.py` (stdlib `urllib`
only, the repository's first network code) downloads a small real, accession-backed *E. coli*
cohort from NCBI (BioProject PRJNA636382) and `test_real_cohort.py`/`make real-cohort-test` runs
the full pipeline against it -- kept out of `make test`/CI since it requires network access,
unlike the rest of this project. Run live against real sequences, it correctly recovered a
persistence signal already known from the source study: both plasmids carried by isolate
`upec_ecpf5` are shared with `upec_ecpf7`, the same recurrent-UTI patient's other episode.

GitHub Actions remains disabled because the available authorization lacks workflow-write scope. The [CI template](ci-template.yml) includes the new tests and can be installed by a maintainer with that permission.
