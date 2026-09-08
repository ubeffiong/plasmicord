#!/usr/bin/env bash
# Stage 0 — make dirs and check dependencies.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config/config.sh"; source "$HERE/lib.sh"
log "Project root: $PROJECT_ROOT"
mkdir -p "$DATA_DIR" "$RESULTS_DIR" "$LOG_DIR" "$TMP_DIR" "$RESULTS_DIR/plasmids"

log "Core deps:"
OK=1
for c in python3 awk; do have "$c" && log "  [ok]   $c" || { warn "  [MISS] $c"; OK=0; }; done

if [[ "$PLASMID_SOURCE" == "mob_recon" ]]; then
    have mob_recon && log "  [ok]   mob_recon" || { warn "  [MISS] mob_recon (PLASMID_SOURCE=mob_recon)"; OK=0; }
else
    log "  [skip] plasmid extraction (PLASMID_SOURCE=$PLASMID_SOURCE)"
fi

if [[ "$DISTANCE_ENGINE" == "mash" ]]; then
    have mash && log "  [ok]   mash" || { warn "  [MISS] mash (DISTANCE_ENGINE=mash)"; OK=0; }
else
    log "  [ok]   kmer engine (pure python, no external dep)"
fi

[[ "$OK" -eq 1 ]] && log "Dependency check PASSED." || warn "Missing deps; see INSTALL.md."
