#!/usr/bin/env bash
# Stage 5 — assemble a plain-text/markdown report of the run.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../config/config.sh"; source "$HERE/lib.sh"
REPORT="$RESULTS_DIR/REPORT.md"
{
  echo "# PlasmiCord: A Chromosome-Aware Plasmid Transmission and Surveillance Framework"
  echo
  echo "Connecting plasmid sharing with chromosomal epidemiology."
  echo
  echo "- Distance engine: \`$DISTANCE_ENGINE\`  |  PU threshold: \`$PU_THRESHOLD\`  |  linkage: \`$LINKAGE\`"
  echo
  if [[ -s "$RESULTS_DIR/plasmid_index.tsv" ]]; then
    echo "- Plasmids indexed: $(($(wc -l < "$RESULTS_DIR/plasmid_index.tsv") - 1))"
  fi
  if [[ -s "$RESULTS_DIR/plasmid_clusters.tsv" ]]; then
    echo "- Plasmid units (PUs): $(($(tail -n +2 "$RESULTS_DIR/plasmid_clusters.tsv" | cut -f2 | sort -u | wc -l)))"
  fi
  echo
  if [[ -s "$RESULTS_DIR/discordance.summary.txt" ]]; then
    echo '## Discordance summary'; echo '```'; cat "$RESULTS_DIR/discordance.summary.txt"; echo '```'; echo
    echo "See \`discordance.crosslinks.tsv\` for the plasmid-mediated links that cross chromosomal clusters."
  fi
  echo
  echo "## Files"
  echo "- \`plasmid_clusters.tsv\` — plasmid_id -> plasmid unit"
  echo "- \`network.edges.tsv\`, \`network.graphml\` — isolate plasmid-sharing network (open GraphML in Cytoscape)"
  echo "- \`network.isolate_units.tsv\` — isolate -> plasmid units carried"
  echo "- \`discordance.crosslinks.tsv\`, \`discordance.summary.txt\` — chromosome-vs-plasmid discordance"
} > "$REPORT"
log "Stage 5 complete -> $REPORT"
cat "$REPORT"
