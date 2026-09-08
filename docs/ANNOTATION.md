# Automatic functional annotation

PlasmiCord's essential profile runs three local stages on each accepted candidate:

| Stage | Engine | Recorded evidence |
|---|---|---|
| Genes | Prokka or Bakta | CDS coordinates, products, cross-references |
| AMR | AMRFinderPlus nucleotide search, with --plus | Reference matches, identity, coverage, method, drug class |
| Mobility | MOB-typer | Candidate-level replicon, relaxase, oriT and mobility typing |

Install the engines and their databases before running. Database downloads are explicit setup operations. Use Linux/WSL for these tools. Activate an environment that supplies their supporting executables, including Perl, BLAST and HMMER. Each configured executable must return a version; missing tools, incompatible outputs and failed commands stop the run with retained logs.

Copy [annotation.example.json](../config/annotation.example.json), replace database paths and release labels, then run:

```sh
plasmicord run --manifest manifest.tsv --metadata metadata.tsv \
  --annotation-config annotation.json --annotation-cache annotation_cache \
  --threads 4 --engine mash --threshold 0.01 --out results_annotated
```

The essential profile requires all three stages. Use `--annotation-profile custom` to explicitly run a subset. Disabled stages remain **not evaluated**. Imported `--features` and automatic annotation are mutually exclusive, preventing accidental mixing of caller results.

Every cache identity contains the candidate sequence, shared-component version, engine version, database content SHA-256, declared database release and configuration. Cached payloads are checksum-verified and original contig IDs are restored. Different labels alone cannot conceal a changed database: actual file contents are fingerprinted once per run. Keep database directories immutable during execution. Cache provenance records the original command; a cache hit need not reproduce raw files in the new output directory.

Outputs include `automatic_features.tsv`, normalized `functional_features.tsv`, `annotation_status.tsv`, `annotation_provenance.json`, `mobility_typing.json`, per-stage logs/raw outputs, and a PlasBench-compatible `plasbench_proteins.tsv`. Report tables distinguish completion with zero matches from an unrun stage.

AMRFinder core-scope calls become headline ARGs only with at least 90% identity and 90% reference coverage by default, and without partial, internal-stop or point-mutation methods. Plus-scope and lower-support calls remain inspectable. This extra reporting rule is not a substitute for AMRFinder's family-specific detection criteria. A nucleotide-only search without an organism setting does not provide comprehensive resistance-mutation assessment. No protein product substring establishes an ARG. MOB typing is aggregate candidate evidence; it does not invent locus coordinates, transfer rates or host-to-host transmission.

Bakta/Prokka gene products do not establish orthology, pathway completeness or resistance phenotype. General CDS products may remain hypothetical. The shared package exports coordinates to PlasBench's zero-based, half-open protein contract while the PlasmiCord feature contract stays one-based inclusive.

Primary tool specifications: [Bakta](https://github.com/oschwengers/bakta), [Prokka](https://github.com/tseemann/prokka), [AMRFinderPlus result interpretation](https://github.com/ncbi/amr/wiki/Interpreting-results), [MOB-suite](https://github.com/phac-nml/mob-suite).
