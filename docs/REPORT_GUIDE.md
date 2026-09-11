# Detailed final-report guide — PlasmiCord v0.4.0

Every completed analysis now produces an offline dashboard, a detailed Markdown companion, a recursive output manifest and a ZIP containing the report with all regular files in its result directory. The report is generated automatically by `plasmicord run` and `plasmicord demo`.

The presentation follows the useful patterns of [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (modular evidence summaries, visible status and offline HTML) and [MultiQC](https://docs.seqera.io/multiqc/reports) (cohort comparison and interactive results). It does not run either tool, use their QC thresholds, or claim to be a MultiQC plugin.

PlasBench's local benchmark report informed the guided findings, explicit unevaluated states, individual-result inspection, accessible status palette and recursive file explorer. PlasmiCord's implementation remains independent: no PlasBench application, templates or runtime are loaded.

## Opening and sharing results

Open **REPORT.html** in a modern browser. It embeds its JavaScript, CSS, chart data, result tables and bounded text previews. It makes no external service requests. It works from `file://` without a web server.

Use **Download complete results ZIP**, then extract the ZIP before opening its REPORT.html. The ZIP includes candidate FASTAs, native annotation outputs and logs that exist in this run directory, tables, provenance, calibration attachments and the final manifest. Keep the extracted directory structure for original-file downloads. A standalone HTML copy still displays its embedded evidence, but cannot supply originals left behind.

Use **Print / save PDF** for a static copy. Printing expands report details and renders every row matching the current table filters, rather than only the current page. Interactive selections and table filters are retained; reset them for an unfiltered print. Graph, chart and gene-track views retain the currently displayed selections and documented limits.

## Dashboard modules

| Module | Detailed contents |
|---|---|
| Overview and findings | Cohort counts, execution/evidence status, review priorities, evidence paths and versioned rule IDs |
| General statistics | Every isolate, accepted/rejected candidate counts, total accepted length, unit count, sharing partners, eligible ARG labels and AMR completion denominator |
| Individual isolate | Full supplied metadata, selected-isolate interpretation, candidate inventory, annotation stages, pair–unit evidence and candidate-specific files |
| Distances and units | Unit sizes, threshold sensitivity, pairwise distance heatmap, exact pair selectors and membership table |
| Sharing and chromosome context | Metadata/ARG filters, chromosome/location/organism colours, optional node-size encoding, clickable nodes/edges, direct-threshold support, a pairwise gene-track comparison for a selected edge, an isolate timeline with distance-labelled edges, cross-cluster observations and filtered-network interpretation |
| Functional cargo | Unit/candidate/category/search selectors, directional contig gene tracks with a circular gene-map toggle, complete feature provenance, observed function prevalence and evaluation coverage |
| Quality and annotation | All candidate decisions, biological evidence, stage completion, versions, cache reuse and raw sensitivity values |
| Calibration | Matching attached calibration's target, training-selected cutoff, current analysis cutoff, holdout confusion counts, metrics and limitations |
| Output explorer | Every regular result file, nested folders, search/module filters, bounded text preview, original-file access/download, type, purpose, byte size, UTC modification time and SHA-256 |
| Run metadata and methods | Analysis and renderer versions, parameters, software/database records, full provenance snapshot, rule explanation and scientific safeguards |

Tables support search, numeric/text sorting, 25-row pagination, keyboard-accessible inspection and filtered TSV export. Exports include all matching rows, not just the displayed page. Spreadsheet formula-like text is escaped in interactive TSV exports; raw output downloads are unchanged.

Charts show quality tiers, candidate size bands, annotation completion, function-category carriers, unit sizes, drug-class carriers, metadata coverage, collection months and threshold sensitivity. Each chart has an SVG download and a TSV export of all its data. Counting rules and denominators are recorded in report_charts.json. Categories and drug classes can overlap; these charts must not be added together as disjoint counts.

## Automated interpretations

Interpretation is **deterministic and evidence-based**, with rule set version 1.0. No hosted model, API key or generated biological narrative is involved.

The cohort rules evaluate:
- Execution completeness and failure state.
- Candidate assignments, actual rejections, low-confidence and uncertain quality.
- Chromosome-typing coverage and observed cross-cluster pairs.
- Sharing observations lacking direct-threshold support.
- AMR observations, completed searches and completion with no eligible calls.
- Gene and mobility caller completion.
- Missing organism, date and location metadata.
- Sensitivity of unit counts to the tested thresholds.
- Attached calibration status, selected cutoff and stated validation target.
- Exact duplicate sequences and annotation-cache reuse.

Every finding records its rule/version, state, interpretation, supporting file paths and suggested review action. **interpretations.json** and **module_summary.tsv** preserve those findings.

Selected-isolate, filtered-network, selected-candidate and exact-pair interpretations update when the corresponding controls change. Cohort-level findings and charts remain explicitly cohort-wide. No filter silently changes the underlying clustering or the analysis threshold.

For example, no supplied chromosome labels produces **not evaluated**, not concordance; a completed AMR search with zero eligible calls is distinguished from an unrun search; single-linkage sharing without a directly qualifying pair is identified for review. A failed run does not convert unassigned candidates into rejected plasmids or fabricate zero completed sharing results.

## Colour semantics

Status is always expressed in text as well as colour:
- **Complete / supported:** execution or evidence coverage is present, not proof of biological accuracy.
- **Review:** missing coverage, questionable quality, sensitivity or evidence requiring examination.
- **Failed:** incomplete execution or explicit rejection where the label states that.
- **Not evaluated:** a required observation or analysis is unavailable.
- **Descriptive:** a measured summary without a biological pass/fail judgment.

An optional blue/orange status palette improves separation. Quality labels retain their exact tier names. Functional categories and network metadata use their own legends. The distance heatmap scales colours to its displayed distance range and is not a quality grade.

## Regenerate an existing result

Refresh the presentation without rerunning reconstruction, annotation, distances or clustering:

```sh
python plasmicord.py report --results results_cohort
```

The report checks recorded analysis-file and report-data hashes before regeneration. Changed/missing analysis evidence is rejected rather than silently mixed with an old report. New extra files can be inventoried, but do not become new biological measurements.

Attach a matching calibration:

```sh
python plasmicord.py report --results results_cohort --calibration calibration
```

The matrix digest, engine and linkage must match; recorded sketch parameters are checked when present. Required calibration JSON, metrics and labels are copied into the result directory. An existing different attachment is not overwritten. The attachment does not change the analysis threshold. Calibration based on reference alignments remains a sequence-relatedness exercise, not transmission validation.

For very large outputs:

```sh
python plasmicord.py report --results results_cohort --no-bundle
```

This regenerates the dashboard without ZIP creation and removes a stale report ZIP if one exists.

## File integrity and preview limits

The explorer includes all regular files within the result directory. Symlinks/junctions leading outside it are not followed. Paths are relative and URL-encoded; text previews are displayed as text, never executed as HTML.

Each text preview is capped at 12,000 bytes, with a total embedded preview budget of 1 MiB. Binary or omitted previews retain original-file download links. The file tree is not truncated by these limits.

The final **output_manifest.json** records file sizes, modification times, purposes and SHA-256 hashes after rendering. The manifest and bundle exclude their own hashes because self-referential hashes cannot be defined. The HTML marks its own final metadata and the changing run-provenance record as available from the final manifest, rather than embedding stale hashes.

**run_provenance.json** contains recursive hashes for stable result files, including native outputs. Its embedded HTML snapshot is explicitly labelled as a rendering-time snapshot. Analysis version, report-renderer version, run completion time and report generation time are separate concepts.

## Display and scientific limits

The heatmap displays at most the first 120 matrix candidates; exact selectors cover that displayed subset. The complete matrix remains downloadable. The network displays at most 200 matching isolates and asks the reader to narrow the search when limited. Bar charts display at most 20 categories while exports retain every category. These are disclosed display limits, not omissions from the underlying evidence.

This release retains PlasmiCord's scientific limits: no inferred transmission direction, no invented chromosome labels or dates, no phenotypic resistance claim, no universal cutoff, no calibrated biological-quality probability, and no pathway completeness inferred from isolated products.

## Verification

The report has been generated and exercised on the seeded demo, the 81-plasmid/32-isolate public-reference cohort with matching calibration, and the two real annotated reference plasmids.

Regression tests cover missing versus completed evidence, actual rejection versus incomplete assignment, failed annotation dashboards, nested outputs, malicious text/filenames, bounded previews, external symlinks, final checksums, complete ZIP contents, calibration matching and regeneration integrity.

```sh
python -m unittest discover -s test -p test_detailed_report.py -v
```

The browser assertions in test/report_browser_checks.js exercise pagination/filtering, selected-isolate navigation, empty network views, exact distances, output previews/download URLs, caller provenance, calibration display and print pagination. They can be evaluated with agent-browser against a generated report; that browser tool is a development QA dependency only.
