#!/usr/bin/env bash
# run_demo.sh -- full framework end-to-end on synthetic sequences, OFFLINE.
# No mob_recon, no Mash: uses the pure-python k-mer distance on tiny plasmids.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
PY="$ROOT/python"
DEMO="$ROOT/results_demo"
rm -rf "$DEMO"; mkdir -p "$DEMO"

echo "### 1) generate synthetic outbreak set"
python3 "$HERE/make_synthetic.py" --outdir "$DEMO"

echo "### 2) all-vs-all distances (pure-python k-mer; Mash used on real data)"
python3 "$PY/kmer_distance.py" --dir "$DEMO/plasmids" --k 15 --out "$DEMO/plasmid_matrix.tsv"

echo "### 3) cluster into plasmid units (single linkage, threshold 0.10)"
python3 "$PY/cluster_plasmids.py" --matrix "$DEMO/plasmid_matrix.tsv" \
    --threshold 0.10 --linkage single --out "$DEMO/plasmid_clusters.tsv"
echo "--- plasmid_clusters.tsv ---"; cat "$DEMO/plasmid_clusters.tsv"

echo "### 4) build isolate-level plasmid-sharing network"
python3 "$PY/build_network.py" --index "$DEMO/plasmid_index.tsv" \
    --clusters "$DEMO/plasmid_clusters.tsv" --metadata "$DEMO/metadata.tsv" \
    --out-prefix "$DEMO/network"
echo "--- network.edges.tsv ---"; cat "$DEMO/network.edges.tsv"

echo "### 5) chromosome-vs-plasmid discordance"
python3 "$PY/discordance.py" --isolate-units "$DEMO/network.isolate_units.tsv" \
    --metadata "$DEMO/metadata.tsv" --out-prefix "$DEMO/discordance"
echo "--- discordance.crosslinks.tsv ---"; cat "$DEMO/discordance.crosslinks.tsv"
echo "--- discordance.summary.txt ---"; cat "$DEMO/discordance.summary.txt"

echo
echo "Demo outputs in: $DEMO"
echo "Open network.graphml in Cytoscape to see the plasmid-sharing network."
