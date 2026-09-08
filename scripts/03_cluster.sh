#!/usr/bin/env bash
# Stage 3 — cluster plasmids into plasmid units (PUs).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config/config.sh"; source "$HERE/lib.sh"; need python3
MATRIX="$RESULTS_DIR/plasmid_matrix.tsv"
[[ -s "$MATRIX" ]] || die "no matrix at $MATRIX; run stage 2."
python3 "$HERE/../python/cluster_plasmids.py" --matrix "$MATRIX" \
    --threshold "$PU_THRESHOLD" --linkage "$LINKAGE" \
    --out "$RESULTS_DIR/plasmid_clusters.tsv"
log "Stage 3 complete -> $RESULTS_DIR/plasmid_clusters.tsv"
