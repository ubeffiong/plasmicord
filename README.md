# PlasmiCord

**PlasmiCord: A Chromosome-Aware Plasmid Transmission and Surveillance Framework**

*Connecting plasmid sharing with chromosomal epidemiology.*

**From reconstruction quality to plasmid-sharing evidence.** This standalone research
framework groups reconstructed plasmids into sequence-defined plasmid units (PUs),
builds isolate-sharing networks, identifies sharing across known chromosomal clusters,
and presents imported functional cargo alongside provenance and quality warnings.

PlasBench asks **“How trustworthy is this reconstruction?”** This project asks
**“Which plasmid units and functions are observed across isolates, places and dates?”**
It does not install, import or call PlasBench. Any reconstruction workflow can supply
the same FASTA-plus-manifest contract. See [project positioning](docs/POSITIONING.md).

**Status: v0.2.0 research prototype.** A sharing link is not proof of direct transmission.
Automated biological annotation and calibrated quality grading remain future work.
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
all use the same normalized manifest and validation gate; they are not automatic
parsers for those tools' native output directories. PlasBench results must be selected
and mapped to this contract just like outputs from MOB-recon, Flye or Unicycler.
Raw FASTQ is rejected; assemble, polish and identify plasmid candidates upstream.

Add `--features functional_features.tsv` for functional cargo reporting. The
[input contract](docs/INPUT_CONTRACT.md) defines coordinates, provenance, quality flags,
headline ARG rules and annotation limitations. Without annotations, the report says
**not evaluated**, rather than claiming genes are absent.

## What the final output looks like

**`REPORT.html` is the primary deliverable:** one self-contained offline report with:

- Overview counts and a prominent dataset/interpretation notice.
- Clickable isolate-sharing network; metadata and ARG filters; node colouring;
  edge details showing distances, shared units, shared eligible ARGs and direct
  threshold support versus links induced through cluster membership.
- Cross-cluster pair–unit evidence table, excluding unknown chromosome clusters.
- **Functional Cargo and AMR Transmission:** unit/isolate selection, category/search
  filters, directional contig gene tracks with zoom, and feature provenance details.
- Unit function prevalence, carrying isolates, organisms, locations and date ranges.
- Quality warnings, duplicate detection, threshold sensitivity, interpretation limits
  and downloadable supporting evidence.

Keep the results folder together for downloads. No network access is required to view
the report. Browser Print can produce a static PDF; retain HTML for interactivity.

| Output | Purpose |
|---|---|
| `REPORT.html`, `REPORT.md` | Interactive final report and plain-text narrative |
| `report_data.json` | Machine-readable report payload |
| `validation.tsv`, `plasmid_index.tsv` | Quality/provenance records and accepted candidates |
| `plasmids/`, `metadata.tsv` | Normalized sequences and metadata used in the run |
| `plasmid_matrix.tsv`, `plasmid_clusters.tsv` | Distances and run-local PU assignments |
| `network.graphml` | Network for Cytoscape or Gephi |
| `network.edges.tsv`, `network.isolate_units.tsv` | Sharing pairs and isolate membership |
| `network.edge_evidence.tsv` | Pair–unit distances, direct support and shared eligible ARGs |
| `discordance.crosslinks.tsv`, `discordance.summary.txt` | Known cross-cluster sharing evidence |
| `functional_features.tsv`, `plasmid_unit_function.tsv` | Normalized features and observed unit functions |
| `threshold_sensitivity.tsv` | PU counts at alternative thresholds |
| `run_provenance.json` | Versions, parameters, input/output checksums and completion/failure status |

## Installation and verification

See [INSTALL.md](INSTALL.md). Core scripts run on Windows, Linux and macOS; MOB-suite
reconstruction and the supplied conda environment target Linux/WSL.

```sh
python -m unittest discover -s test -p test_standalone.py -v
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
or failed extraction as evidence that an isolate has no plasmids.

The earlier `transmission` command and `transmission.py` launcher remain compatible.

## License

MIT. See [LICENSE](LICENSE). Project scope and methods are described in
[POSITIONING.md](docs/POSITIONING.md) and [METHODS.md](docs/METHODS.md).
