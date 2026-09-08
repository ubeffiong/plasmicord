#!/usr/bin/env bash
# run_all.sh — run the whole framework. Usage: run_all.sh [stage ...]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config/config.sh"; source "$HERE/lib.sh"
STAGES=("$@"); [[ ${#STAGES[@]} -eq 0 ]] && STAGES=(0 1 2 3 4 5)
declare -A S=( [0]=00_setup.sh [1]=01_extract_plasmids.sh [2]=02_distances.sh
               [3]=03_cluster.sh [4]=04_network.sh [5]=05_report.sh )
log "############ PlasmiCord ############"
log "stages: ${STAGES[*]}"
for s in "${STAGES[@]}"; do
    scr="${S[$s]:-}"; [[ -z "$scr" ]] && die "unknown stage '$s'"
    log ">>>>> STAGE $s : $scr"; bash "$HERE/$scr"
done
log "############ DONE ############"
