#!/usr/bin/env bash
# =============================================================================
# Plasmid-aware transmission framework — configuration
# Edit this, then run scripts/run_all.sh
# =============================================================================
export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

export DATA_DIR="$PROJECT_ROOT/data"
export RESULTS_DIR="$PROJECT_ROOT/results"
export LOG_DIR="$PROJECT_ROOT/logs"
export TMP_DIR="${TMP_DIR:-$PROJECT_ROOT/tmp}"

# --- Input sample sheet ------------------------------------------------------
# TSV with header. REQUIRED columns:
#   isolate_id      unique short id (no spaces, no '__')
#   assembly_path   path to that isolate's genome assembly FASTA
# OPTIONAL columns (used if present):
#   chromosomal_cluster   cgMLST/SNP cluster id -> enables discordance analysis
#   date, location        annotate network nodes
export METADATA="$PROJECT_ROOT/config/metadata.tsv"

# --- Compute -----------------------------------------------------------------
export THREADS="${THREADS:-4}"

# --- Plasmid extraction ------------------------------------------------------
# How to get per-isolate plasmid sequences from each assembly:
#   mob_recon  : run MOB-suite mob_recon on each assembly (default, recommended)
#   precomputed: assemblies ARE already single-plasmid FASTAs, or you dropped
#                per-plasmid FASTAs into results/plasmids yourself -> skip stage 1
export PLASMID_SOURCE="${PLASMID_SOURCE:-mob_recon}"

# --- Distance engine ---------------------------------------------------------
#   mash  : sketch + all-vs-all (fast, recommended for real datasets)
#   kmer  : pure-python exact k-mer distance (no deps; SMALL datasets/demo only)
export DISTANCE_ENGINE="${DISTANCE_ENGINE:-mash}"
export MASH_SKETCH_SIZE="${MASH_SKETCH_SIZE:-10000}"       # mash -s ; larger = more accurate, slower
export MASH_KMER="${MASH_KMER:-21}"                 # mash -k
export KMER_K="${KMER_K:-21}"                    # used when DISTANCE_ENGINE=kmer

# --- Clustering --------------------------------------------------------------
# Threshold on distance to call two plasmids the same "plasmid unit" (PU).
# For Mash distance, ~0.05 ≈ 95% ANI. Tune for your taxon; report what you used.
export PU_THRESHOLD="${PU_THRESHOLD:-0.05}"
export LINKAGE="${LINKAGE:-single}"             # single | complete

if [[ ! -f "$METADATA" ]]; then
    echo "WARNING: metadata sheet not found at $METADATA" >&2
fi
