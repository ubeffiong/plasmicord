# Input & output formats

## Input: config/metadata.tsv (tab-separated, has header)
Required: `isolate_id`, `assembly_path`.
Optional: `chromosomal_cluster` (enables discordance), `date`, `location`.
Rules: `isolate_id` has no spaces and no `__` (double underscore is reserved as the
isolate/plasmid separator). Lines starting with `#` and blank lines are ignored.

Example:
```
isolate_id	assembly_path	chromosomal_cluster	date	location
iso001	/data/asm/iso001.fasta	CC12	2025-01-14	HospitalA
iso002	/data/asm/iso002.fasta	CC12	2025-01-20	HospitalA
iso003	/data/asm/iso003.fasta	CC44	2025-02-02	HospitalB
```

## Intermediate: results/plasmid_index.tsv
`plasmid_id <TAB> isolate_id <TAB> length`. Built by stage 1. If you use
`PLASMID_SOURCE=precomputed`, you create this yourself and drop matching
`<plasmid_id>.fasta` files (one record each) into `results/plasmids/`.

## Intermediate: results/plasmid_matrix.tsv
Square symmetric distance matrix. Header row: `item` then all ids; each row: id then
distances. Diagonal ~0. Produced by Mash (via `mash_to_matrix.py`) or `kmer_distance.py`.

## Outputs
- `plasmid_clusters.tsv` — `plasmid_id <TAB> plasmid_unit`
- `network.isolate_units.tsv` — `isolate_id <TAB> n_plasmid_units <TAB> comma,list,of,PUs`
- `network.edges.tsv` — `source <TAB> target <TAB> weight <TAB> shared_units`
- `network.graphml` — undirected network with node/edge attributes
- `discordance.crosslinks.tsv` — `plasmid_unit, isolate_a, chrom_cluster_a, isolate_b, chrom_cluster_b`
- `discordance.summary.txt`, `REPORT.md` — human-readable summaries
