#!/usr/bin/env bash
# Stage 2 — all-vs-all plasmid distances -> results/plasmid_matrix.tsv
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config/config.sh"; source "$HERE/lib.sh"
PDIR="$RESULTS_DIR/plasmids"
MATRIX="$RESULTS_DIR/plasmid_matrix.tsv"

shopt -s nullglob
fastas=("$PDIR"/*.fasta)
shopt -u nullglob
[[ ${#fastas[@]} -ge 1 ]] || die "need >=1 plasmid FASTA in $PDIR (found ${#fastas[@]}). Run stage 1."

if [[ "$DISTANCE_ENGINE" == "mash" ]]; then
    need mash
    log "Mash sketching (${#fastas[@]} plasmids, k=$MASH_KMER, s=$MASH_SKETCH_SIZE) ..."
    mash sketch -k "$MASH_KMER" -s "$MASH_SKETCH_SIZE" -o "$RESULTS_DIR/plasmids_sketch" \
        "${fastas[@]}" > "$LOG_DIR/mash_sketch.log" 2>&1 || die "mash sketch failed"
    log "Mash all-vs-all dist ..."
    mash dist "$RESULTS_DIR/plasmids_sketch.msh" "$RESULTS_DIR/plasmids_sketch.msh" \
        > "$RESULTS_DIR/mash_dist.tsv" 2> "$LOG_DIR/mash_dist.log" || die "mash dist failed"
    python3 "$HERE/../python/mash_to_matrix.py" --dist "$RESULTS_DIR/mash_dist.tsv" --out "$MATRIX"
else
    log "pure-python k-mer distances (k=$KMER_K) ..."
    python3 "$HERE/../python/kmer_distance.py" --dir "$PDIR" --k "$KMER_K" --out "$MATRIX"
fi
log "Stage 2 complete -> $MATRIX"
