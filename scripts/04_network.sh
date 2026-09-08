#!/usr/bin/env bash
# Stage 4 — isolate-level plasmid-sharing network + discordance vs chromosome.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config/config.sh"; source "$HERE/lib.sh"; need python3
INDEX="$RESULTS_DIR/plasmid_index.tsv"
CLUST="$RESULTS_DIR/plasmid_clusters.tsv"
[[ -s "$INDEX" && -s "$CLUST" ]] || die "need $INDEX and $CLUST; run stages 1-3."

python3 "$HERE/../python/build_network.py" --index "$INDEX" --clusters "$CLUST" \
    --metadata "$METADATA" --out-prefix "$RESULTS_DIR/network"

# discordance only if chromosomal_cluster column present
if head -1 "$METADATA" | tr '\t' '\n' | grep -qx "chromosomal_cluster"; then
    python3 "$HERE/../python/discordance.py" \
        --isolate-units "$RESULTS_DIR/network.isolate_units.tsv" \
        --metadata "$METADATA" --out-prefix "$RESULTS_DIR/discordance"
else
    warn "no 'chromosomal_cluster' column in metadata; skipping discordance analysis."
fi
log "Stage 4 complete."
