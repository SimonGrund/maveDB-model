#!/usr/bin/env bash
set -x

shopt -s nullglob

# Base folders
INPUT_ROOT="/app/Shared/Input"
OUTPUT_ROOT="/app/Shared/Output"

# Fixed file names inside each subfolder
CLINICAL_FILE_NAME="clinical.csv"
DEPMAP_FILE_NAME="depmap.csv"

# Infinite loop
while true; do

  echo "Loop running..."

  for dir in "$INPUT_ROOT"/*; do
    # Skip if not a directory
    [ -d "$dir" ] || continue

    SUBDIR_NAME=$(basename "$dir")

    CLINICAL_FILE="$dir/$CLINICAL_FILE_NAME"
    DEPMAP_FILE="$dir/$DEPMAP_FILE_NAME"
    OUTPUT_DIR="$OUTPUT_ROOT/$SUBDIR_NAME"

    # Check that both input files exist
    if [ ! -f "$CLINICAL_FILE" ]; then
      echo "Skipping $SUBDIR_NAME – missing $CLINICAL_FILE_NAME"
      continue
    fi

    if [ ! -f "$DEPMAP_FILE" ]; then
      echo "Skipping $SUBDIR_NAME – missing $DEPMAP_FILE_NAME"
      continue
    fi

    echo "Processing $SUBDIR_NAME..."
    Rscript Code/run_pipeline.R "$CLINICAL_FILE" "$DEPMAP_FILE" "$OUTPUT_DIR"

	rm -f "$dir"/*
	rmdir "$dir"

  done

  # Wait before next scan
  sleep 5
done

echo "Loop ended"


