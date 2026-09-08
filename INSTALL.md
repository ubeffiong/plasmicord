# PlasmiCord installation

## Portable standalone core

Python 3.10 or later is sufficient on Windows, Linux or macOS:

```sh
python plasmicord.py demo --out results_demo
```

Optionally install the CLI:

```sh
python -m pip install .
plasmicord check
plasmicord demo --out results_demo_installed
```

Open the generated `results_demo/results/REPORT.html`. The browser needs no server,
internet connection or plugins. Keep the result folder for TSV/JSON/GraphML downloads.

## Mash and optional reconstruction tools (Linux/WSL)

With conda or mamba installed:

```sh
conda env create -f env/environment.yml
conda activate plasmicord
python -m pip install .
plasmicord check
plasmicord demo --engine mash --out results_mash_demo
```

The environment includes Mash, MOB-suite and Python. MOB-suite and its databases are
only needed to reconstruct candidates from assemblies, not for imported FASTA analysis.
Database installation can be substantial; initialize and version the database according
to the installed MOB-suite release before reconstruction. The core does not download
databases during analysis. An explicit configuration runs the installed annotators;
see [automatic annotation](docs/ANNOTATION.md).

```sh
plasmicord run --manifest manifest.tsv --metadata metadata.tsv \
  --engine mash --threshold 0.01 --out results_cohort
```

Use an appropriately calibrated threshold; the example is not a recommended universal
value. Use `--engine kmer` explicitly for small inputs if Mash is unavailable.

## Troubleshooting

- Missing Mash: activate the environment containing `mash`, then run `plasmicord check`.
- Existing output directory: choose a fresh directory. Completed or failed runs are
  preserved so an interrupted run cannot silently mix with an older result.
- Unknown chromosome clusters: sharing networks remain available; those isolates are
  excluded from cross-cluster counts.
- Missing functional rows: configure `--annotation-config` or import `--features`.
  Inspect stage completion before interpreting zero hits.
- Rejected short sequences: inspect `validation.tsv` and the study's `--min-length`.
- Broken links after sharing HTML: distribute the whole results folder for downloads.
- Legacy extraction failure: inspect the MOB log and resolve it. The shell pipeline
  now stops rather than treating failed isolates as having no plasmids.

The older `bash scripts/run_all.sh` and `bash test/run_demo.sh` workflows remain for
compatibility and produce legacy outputs; see the README before using them.
