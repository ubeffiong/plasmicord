# Native imports

All importers write independent FASTAs, a manifest, isolate metadata, validation records and `import_provenance.json`. They preserve source versions and input checksums, then use the same manifest gate. Run `plasmicord run` on the generated files.

```sh
plasmicord import mob-recon --input mob_output --isolate-id ISO001 --tool-version 3.1.9 --out imported_mob
plasmicord import flye --input flye_output --isolate-id ISO001 --record-ids contig_2,contig_5 --out imported_flye
plasmicord import unicycler --input unicycler_output --isolate-id ISO001 --record-ids 2,5 --out imported_hybrid
plasmicord import plasbench --input selected_candidate/ISO001 --isolate-id ISO001 --out imported_plasbench
plasmicord import generic-fasta --input candidate.fasta --isolate-id ISO001 --group-contigs --out imported_fasta
plasmicord import tadrep --input tadrep_output/ISO001 --isolate-id ISO001 --tool-version 1.0.0 --out imported_tadrep
plasmicord run --manifest imported_mob/manifest.tsv --metadata imported_mob/metadata.tsv --threshold 0.01 --out results_mob
```

| Adapter | Required native inputs / interpretation |
|---|---|
| MOB-recon | `contig_report.txt` and `plasmid_*.fasta`; preserves reconstructed multi-contig groups and checks their plasmid classification |
| Flye | `assembly.fasta` and `assembly_info.txt`; explicit plasmid record selection required; preserves coverage, graph and circularity evidence |
| Unicycler | `assembly.fasta`; explicit record selection required; preserves native circular headers and graph provenance |
| PlasBench | One selection report and selected plasmid FASTA; ambiguous multi-record candidates require explicit selection or grouping |
| Generic FASTA | Single candidate, explicit selected records, or deliberate `--group-contigs` for one fragmented candidate |
| [TaDReP](https://github.com/oschwengers/tadrep) | One `<sample>-summary.tsv` and one or more `<sample>-<reference>-pseudo.fna`; each reconstructed reference plasmid is one candidate automatically; preserves coverage/identity as the conservative floor and alignment length as the sum across contributing contig rows |

For batches, replace `--input` and `--isolate-id` with `--samples samples.tsv`. Required columns: `isolate_id,input_path` (tab-separated); optional `record_ids,source_tool_version,chromosomal_cluster,date,location,organism`. Paths resolve from the samples file. Record IDs are comma-separated.

Do not select assembly contigs solely because they are circular: chromosomes can also be circular. These adapters consume completed native outputs; they do not assemble raw FASTQ. PlasBench is an optional file format, with no import or runtime dependence on that application. Its recovery scores are source context, not transmission confidence.

Native specifications: [MOB-suite](https://github.com/phac-nml/mob-suite), [Flye](https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md), [Unicycler](https://github.com/rrwick/Unicycler), [TaDReP](https://github.com/oschwengers/tadrep).

Importing `mob-recon` output does **not** automatically populate MOB-suite's own
`primary_cluster_id`/`secondary_cluster_id` -- those come only from `mob_cluster`
(a separate MOB-suite stage this adapter does not run) and can be attached afterward via
`--external-typing` on `plasmicord run`; see [input contract](INPUT_CONTRACT.md).

TaDReP reports no circularity/completeness data at all; the `tadrep` adapter always records
`circularity_status=unresolved`, never inferring it. TaDReP's own `plasmids.info.tsv` (cohort
plasmid characterization, not imported here) has a confirmed upstream key-mismatch bug in its
`CDS` column (`characterize.py` sets `plasmid['cds']`, but `detect.py`'s cohort-info writer
reads `['cdrs']`) -- do not rely on that column if wiring it in later.
