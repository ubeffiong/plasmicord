#!/usr/bin/env bash
# Stage 1 — extract per-isolate plasmid FASTAs into results/plasmids and build
# results/plasmid_index.tsv (plasmid_id, isolate_id, length).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config/config.sh"; source "$HERE/lib.sh"
ADAPT="$HERE/../adapters"
PDIR="$RESULTS_DIR/plasmids"; mkdir -p "$PDIR"
INDEX="$RESULTS_DIR/plasmid_index.tsv"

[[ -f "$METADATA" ]] || die "metadata sheet not found: $METADATA"

if [[ "$PLASMID_SOURCE" == "precomputed" ]]; then
    log "PLASMID_SOURCE=precomputed: expecting per-plasmid FASTAs already in $PDIR"
    [[ -s "$INDEX" ]] || die "no $INDEX found; create it (plasmid_id<TAB>isolate_id<TAB>length) or use mob_recon."
    exit 0
fi

need mob_recon
# Avoid silently mixing a previous cohort with a fresh extraction.
[[ -z "$(find "$PDIR" -maxdepth 1 -type f -name '*.fasta' -print -quit)" ]] || die "plasmids directory is not empty; use a fresh RESULTS_DIR"
# fresh index
printf "plasmid_id\tisolate_id\tlength\n" > "$INDEX"

while IFS=$'\t' read -r ISO ASM; do
    [[ -z "${ISO:-}" ]] && continue
    [[ -s "$ASM" ]] || die "assembly missing for $ISO: $ASM"
    log "=== extract plasmids: $ISO ==="
    bash "$ADAPT/extract_mob_recon.sh" "$ISO" "$ASM" "$PDIR" "$INDEX" "$THREADS" \
        "$LOG_DIR/${ISO}.mob_recon.log" || die "extraction failed for $ISO; cohort is incomplete"
done < <(read_isolates "$METADATA")

n=$(($(wc -l < "$INDEX") - 1))
log "Stage 1 complete: $n plasmid(s) indexed -> $INDEX"
