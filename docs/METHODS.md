# Methods and interpretation boundaries

## Inputs and validation

The core accepts reconstructed candidate FASTAs and isolate metadata through the
[input contract](INPUT_CONTRACT.md). Source modes share technical checks and provenance
normalization. A technically valid candidate is not necessarily biologically a plasmid.
Source-provided circularity and quality claims remain distinguishable from evaluation.
Reconstruction is upstream; the optional legacy shell workflow uses MOB-recon.

## Distances

Mash sketches each candidate and computes all-versus-all distances. The fallback
computes exact canonical k-mer sets and transforms Jaccard similarity J using
`D = -ln(2J/(1+J))/k`, with distance 1 for no shared k-mers. Ambiguous k-mers are
excluded, and contig boundaries are preserved. Candidate order and normalized
sequence checksums are retained. The converter rejects missing pairs instead of
inventing zero distances. Matrices must be square, symmetric, finite, bounded in
[0,1], uniquely labelled and zero on the diagonal.

Mash distance is an approximate mutation-distance model, not direct alignment identity,
and size differences/accessory regions affect similarity. See the
[official distance formulation](https://mash.readthedocs.io/en/latest/distances.html)
and [Mash methods paper](https://doi.org/10.1186/s13059-016-0997-x).

## Plasmid units

Single linkage takes connected components of the graph with distances <= threshold.
Complete linkage agglomerates the pair of clusters with the smallest maximum
between-cluster distance, provided it remains <= threshold. Lexicographic member
names break equal-distance ties deterministically. The new CLI defaults to complete
linkage; the legacy scripts and synthetic demo use single linkage.

PU IDs are assigned by sorted smallest member names. They are reproducible for a
fixed dataset, but **run-local**, not stable registry identifiers across changing cohorts.
The CLI requires a threshold and reports sensitivity at zero, half, selected and
double threshold (bounded by 1). This describes sensitivity, not threshold calibration.

## Network and chromosome discordance

Each metadata isolate is a node; edges connect isolates carrying a common PU.
Weight counts distinct shared units. Lack of an indexed candidate is not proof of
plasmid absence. Edge evidence records the minimum candidate distance within each
shared unit and whether that pair directly meets the threshold. Single-linkage
membership can otherwise induce links through intermediate members.

Discordance requires two known, comparable chromosome-cluster labels that differ.
Unknown labels are excluded. Outputs distinguish unique isolate pairs from
pair–unit links. A zero count is no detected discordance among typed isolates, not
proof of concordance or no transmission. Dates are displayed as observations;
the graph remains undirected. Cluster labels across organisms/schemes require
upstream namespace normalization.

## Functional cargo

The core validates and imports normalized features; it does not run gene calling,
specialized AMR detection or orthology assignment. Gene tracks show supplied coordinates
and strand. Headline rules exclude loose/low-confidence/incomplete-provenance ARG calls.
Functional prevalence deduplicates candidate carriers per label/category and uses all
accepted candidates in the unit as denominator. It does not establish orthology,
biological gene absence, validated core/accessory genes or complete metabolic modules.

## Reproducibility and limits

Run provenance contains framework/Python/Mash versions where relevant, parameters,
input checksums, sequence checksums and output hashes. A fresh output folder is
required; failures leave a failed status. HTML is offline and JSON/TSV/GraphML preserve
underlying evidence. No PlasBench installation is required.

Distance storage is quadratic; the exact fallback is capped at 200 candidates.
Complete linkage and dense isolate-pair enumeration can be expensive. Scalable
candidate screening, epidemiological inference, alignment-based confirmation,
biological QC, automatic annotation and real-cohort calibration remain open work.
