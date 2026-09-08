#!/usr/bin/env bash
# Adapter: run mob_recon on ONE assembly and emit each reconstructed plasmid as a
# separate FASTA named <isolate>__<clusterid>.fasta into the shared plasmids dir.
# Also appends rows to a per-isolate index passed in.
#
# Usage: extract_mob_recon.sh <isolate_id> <assembly_fasta> <plasmids_out_dir> <index_tsv> <threads> <logfile>
set -euo pipefail
ISO="$1"; ASM="$2"; PDIR="$3"; INDEX="$4"; THREADS="$5"; LOG="$6"

command -v mob_recon >/dev/null 2>&1 || { echo "mob_recon not found" >&2; exit 1; }
[[ -s "$ASM" ]] || { echo "assembly not found: $ASM" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

mob_recon --infile "$ASM" --outdir "$WORK/out" --num_threads "$THREADS" --force \
    > "$LOG" 2>&1 || { echo "mob_recon failed for $ISO (see $LOG)" >&2; exit 1; }

shopt -s nullglob
n=0
for f in "$WORK"/out/plasmid_*.fasta; do
    # cluster id from filename: plasmid_<clusterid>.fasta
    base="$(basename "$f")"; cid="${base#plasmid_}"; cid="${cid%.fasta}"
    pid="${ISO}__${cid}"
    dest="$PDIR/${pid}.fasta"
    # Preserve contig boundaries; concatenating creates unsupported sequence joins.
    awk -v h=">$pid" 'BEGIN{n=0} /^>/{n++; print h "_contig" n; next} {print}' \
        "$f" > "$dest"
    # length
    len=$(grep -v '^>' "$dest" | tr -d '\n' | wc -c | tr -d ' ')
    printf "%s\t%s\t%s\n" "$pid" "$ISO" "$len" >> "$INDEX"
    n=$((n+1))
done
shopt -u nullglob

echo "[extract_mob_recon] $ISO: $n plasmid(s)" >&2
exit 0
