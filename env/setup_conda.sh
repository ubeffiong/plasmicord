#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if command -v mamba >/dev/null 2>&1; then SOLVER=mamba
elif command -v conda >/dev/null 2>&1; then SOLVER=conda
else echo "ERROR: install Miniforge first (see INSTALL.md)." >&2; exit 1; fi
echo "[setup_conda] using $SOLVER"
"$SOLVER" env create -f "$HERE/environment.yml" || "$SOLVER" env update -f "$HERE/environment.yml"
echo "[setup_conda] done -> conda activate plasmicord"
