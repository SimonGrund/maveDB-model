#!/usr/bin/env bash
set -euo pipefail

# Entrypoint for container and local runs
# Usage:
#   ./run.sh [clin_data_path] [depmap_export_path]
# Defaults:
#   clin_data_path: Data/clin_data.csv
#   depmap_export_path: Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

CLIN_PATH=${1:-"Data/clin_data.csv"}
DEPMAP_PATH=${2:-"Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv"}

echo "Running pipeline with:"
echo "  clin_data:   $CLIN_PATH"
echo "  depmap data: $DEPMAP_PATH"

exec Rscript Code/run_pipeline.R "$CLIN_PATH" "$DEPMAP_PATH"
