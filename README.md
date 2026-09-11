# PlasmiCord

**PlasmiCord: A Chromosome-Aware Plasmid Transmission and Surveillance Framework**

*Connecting plasmid sharing with chromosomal epidemiology.*

**From reconstruction quality to plasmid-sharing evidence.** This standalone research
framework groups reconstructed plasmids into sequence-defined plasmid units (PUs),
builds isolate-sharing networks, identifies sharing across known chromosomal clusters,
and presents automatically detected or imported functional cargo alongside provenance and quality evidence.

PlasBench asks **“How trustworthy is this reconstruction?”** This project asks
**“Which plasmid units and functions are observed across isolates, places and dates?”**
It does not install, import or call PlasBench. Any reconstruction workflow can supply
the same FASTA-plus-manifest contract. See [project positioning](docs/POSITIONING.md).

**Status: v0.4.0 research prototype.** A sharing link is not proof of direct transmission.
Automatic annotation, biological evidence grading, native imports and study-specific
calibration are implemented. Quality tiers remain heuristic research rules.
See the [review and gap assessment](docs/REVIEW.md) for implemented versus planned features.

## Run immediately

Python 3.10+ is sufficient for the small synthetic demo; no third-party Python packages,
bioinformatics databases, browser server, or PlasBench installation are required.

```sh
git clone https://github.com/ubeffiong/plasmicord.git
cd plasmicord
python plasmicord.py demo --out results_demo
```

Open **`results_demo/results/REPORT.html`** in a browser. The demo uses seeded random
sequences and clearly labelled illustrative gene annotations, not biological detections.
Expected: 5 isolates, 5 plasmids, 2 PUs, 4 sharing pairs, 3 cross-cluster pairs.

Optional installation provides the `plasmicord` command from any directory:

```sh
python -m pip install .
plasmicord check
plasmicord demo --out results_demo_installed
```

The exact k-mer engine is intended for small sets (maximum 200 plasmids), not a
replacement for scalable sketching. Real-data runs default to Mash and require an
explicit study-selected threshold. A threshold is not a transmission probability.

## Run your reconstructed plasmids

Provide a tab-separated manifest with three required columns:

```text
isolate_id	plasmid_id	fasta_path
ISO001	ISO001_p1	plasmids/ISO001_p1.fasta
ISO002	ISO002_p1	plasmids/ISO002_p1.fasta
```

Paths are relative to the manifest. Each FASTA represents one plasmid candidate;
multiple contigs in that file represent a fragmented candidate, not multiple plasmids.
Keep separate candidate plasmids in separate files.

Provide a separate metadata TSV containing every isolate, including those without
indexed plasmids. `isolate_id` is required; `chromosomal_cluster`, ISO-format `date`,
`location`, and `organism` are optional. Cluster labels must use a common, documented
typing scheme and be namespaced where organisms/schemes differ.

```sh
plasmicord run --manifest manifest.tsv --metadata metadata.tsv \
  --mode precomputed --engine mash --threshold 0.01 --linkage complete \
  --out results_cohort
```

`0.01` is an example, **not a validated universal cutoff**. Evaluate sensitivity and
confirm prioritized links using sequence alignment, assembly evidence and epidemiology.

Use `--mode longread`, `hybrid`, or `plasbench` to describe the source. These modes
all use the same normalized manifest and validation gate. `plasmicord import` adapts
native MOB-recon, Flye, Unicycler, PlasBench, TaDReP and generic FASTA outputs;
see [native imports](docs/IMPORTS.md).
Raw FASTQ is rejected; assemble, polish and identify plasmid candidates upstream.

Add `--features functional_features.tsv` for functional cargo reporting. The
[input contract](docs/INPUT_CONTRACT.md) defines coordinates, provenance, quality flags,
headline ARG rules and annotation limitations. Without annotations, the report says
**not evaluated**, rather than claiming genes are absent.

For automatic annotation, add `--annotation-config annotation.json` with configured
Prokka/Bakta, AMRFinderPlus and MOB-typer databases. See [annotation and caching](docs/ANNOTATION.md),
[biological quality](docs/QUALITY.md), and [the reusable component](shared_annotation/README.md).

`plasmicord calibrate` fits training labels and evaluates a separate holdout, rejecting
group/isolate/sequence leakage. See [calibration](docs/CALIBRATION.md) and
[public-reference validation](docs/validation/README.md). Alignment agreement does not establish transmission.

Add `--external-typing external_typing.tsv` to cross-reference opaque external identifiers
(MOB-suite cluster IDs, COPLA PTU, predicted host range, PlasmidFinder Inc-types, or a
sequence-based transmissibility classifier's score/call) without letting them influence
PlasmiCord's own quality tiers. `plasmicord population-summary --results results_cohort`
aggregates a completed run into cohort-level plasmid-unit and metadata-dimension summary
tables; see [input contract](docs/INPUT_CONTRACT.md) and [population summary](docs/POPULATION_SUMMARY.md).
Add `--multilayer-network` to `plasmicord run` for an additional per-plasmid-unit GraphML/TSV
export, each edge tagged with whether it stays within or crosses a known chromosomal cluster.
Add `--size-correction-per-percent` to loosen the clustering threshold for plasmid pairs with
large length differences (disabled by default; see [input contract](docs/INPUT_CONTRACT.md)).

## What the final output looks like

**`REPORT.html` is the primary deliverable:** a detailed offline dashboard with a fixed
module sidebar, automated evidence-linked interpretations, colour-coded status,
cohort charts, searchable/sortable/paginated tables and individual-isolate reports.
See the [detailed report guide](docs/REPORT_GUIDE.md) and
[how to read the output](docs/INTERPRETATION.md).

It includes:

- Overview counts and a prominent dataset/interpretation notice.
- Clickable isolate-sharing network; metadata and ARG filters; node colouring and
  optional node-size encoding (sharing degree or candidate-unit count); edge details
  showing distances, shared units, shared eligible ARGs and direct threshold support
  versus links induced through cluster membership; a pairwise gene-track comparison
  for the two plasmids behind a selected edge (not sequence alignment).
- An isolate timeline plotting collection dates with distance-labelled sharing edges;
  isolates without a usable date are listed separately, not omitted.
- Cross-cluster pair–unit evidence table, excluding unknown chromosome clusters.
- **Functional Cargo and AMR Transmission:** unit/isolate selection, category/search
  filters, directional contig gene tracks with zoom (with a circular gene-map toggle,
  one ring per contig), and feature provenance details.
- Unit function prevalence, carrying isolates, organisms, locations and date ranges.
- Quality warnings, duplicate detection, annotation completion and threshold sensitivity.
- Charts for quality, lengths, annotation coverage, functional/drug-class carriers,
  metadata, collection months and unit sizes (with a log/linear scale toggle);
  exact pairwise-distance heatmap.
- Dynamic interpretations for the selected isolate, candidate, pair and filtered network.
- A searchable recursive file tree with metadata, checksums, safe text previews and
  original downloads, including nested annotation outputs and logs.
- A complete `REPORT_BUNDLE.zip`, final file manifest, chart data and detailed Markdown report.

Regenerate the report without repeating analysis, optionally attaching matching calibration:

```sh
plasmicord report --results results_cohort --calibration calibration
```

The command verifies recorded result hashes. Calibration must match the matrix,
engine and linkage; attaching it does not change the analysis threshold. Failed
in-run stages also generate diagnostic HTML when the available evidence can be rendered.

Keep the results folder together for downloads. No network access is required to view
the report. Browser Print can produce a static PDF; retain HTML for interactivity.

| Output | Purpose |
|---|---|
| `REPORT.html`, `REPORT.md` | Detailed interactive dashboard and full plain-text narrative |
| `REPORT_BUNDLE.zip`, `output_manifest.json` | Portable result bundle and recursive final inventory/checksums |
| `interpretations.json`, `module_summary.tsv` | Versioned rules, evidence-linked findings and review actions |
| `sample_summary.tsv`, `report_charts.json` | Per-isolate statistics and exact chart values/denominators |
| `report_data.json` | Machine-readable report payload |
| `validation.tsv`, `plasmid_index.tsv` | Quality/provenance records and accepted candidates |
| `biological_quality.tsv` | Identity, sequence, marker, graph and supplied read/contamination evidence |
| `annotation_status.tsv`, `annotation_provenance.json` | Per-candidate stage completion, versions and cache provenance |
| `mobility_typing.json`, `plasbench_proteins.tsv` | Aggregate typing and reusable protein-coordinate export |
| `plasmids/`, `metadata.tsv` | Normalized sequences and metadata used in the run |
| `plasmid_matrix.tsv`, `plasmid_clusters.tsv` | Distances and run-local PU assignments, incl. per-unit threshold margin |
| `network.graphml` | Network for Cytoscape or Gephi |
| `network.edges.tsv`, `network.isolate_units.tsv` | Sharing pairs and isolate membership |
| `network.edge_evidence.tsv` | Pair–unit distances, direct support, threshold margin and shared eligible ARGs |
| `discordance.crosslinks.tsv`, `discordance.summary.txt` | Known cross-cluster sharing evidence |
| `functional_features.tsv`, `plasmid_unit_function.tsv` | Normalized features and observed unit functions |
| `threshold_sensitivity.tsv` | PU counts at alternative thresholds |
| `typing_crossreference.tsv` | Present only with `--external-typing`; opaque external identifiers (MOB-suite cluster IDs, COPLA PTU, host range), never used to compute quality tiers |
| `containment_candidates.tsv` | Length/similarity heuristic pairs; not alignment-confirmed containment |
| `population_summary.pu_level.tsv`, `population_summary.metadata_dimension.tsv` | Optional cohort-level aggregation from `plasmicord population-summary` |
| `network.multilayer.graphml`, `network.multilayer_edges.tsv` | Present only with `--multilayer-network`; one edge per (isolate pair, plasmid unit), tagged with `cluster_relation` |
| `run_provenance.json` | Versions, parameters, input/output checksums and completion/failure status |

## Installation and verification

See [INSTALL.md](INSTALL.md). Core scripts run on Windows, Linux and macOS; MOB-suite
reconstruction and the supplied conda environment target Linux/WSL.

```sh
python -m unittest discover -s test -p test_standalone.py -v
python -m unittest discover -s test -p test_gap_workflows.py -v
python -m unittest discover -s test -p test_detailed_report.py -v
python -m unittest discover -s test -p test_documentation_links.py -v
python test/test_clustering.py
python test/test_network_discordance.py
python plasmicord.py demo --out results_check
```

Tests establish software behaviour on fixtures; they do not establish biological
accuracy. The [review](docs/REVIEW.md) records validation performed for this release.
An optional Windows/Linux [CI template](docs/ci-template.yml) is included. GitHub
Actions is not enabled; installing the template requires workflow write permission.

## Legacy assembly workflow

The original Ubuntu shell stages remain in `scripts/` and `adapters/` for compatibility:
`bash scripts/run_all.sh` runs MOB-recon, Mash, clustering and a basic Markdown report.
This older route **does not use the new common validation or HTML/functional report**.
Prefer producing a manifest from reconstructed candidate FASTAs and using `plasmicord run`.
The legacy `config/metadata.tsv` contains placeholders only. Do not treat an empty run
or failed extraction as evidence that an isolate has no plasmids. See
[legacy input/output formats](docs/INPUT_FORMATS.md) for this older pipeline's file layout.

The earlier `transmission` command and `transmission.py` launcher remain compatible.

## License

MIT. See [LICENSE](LICENSE). Project scope and methods are described in
[POSITIONING.md](docs/POSITIONING.md) and [METHODS.md](docs/METHODS.md).
