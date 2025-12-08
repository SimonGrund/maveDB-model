#!/usr/bin/env bash
set -x

shopt -s nullglob

# Base folders
INPUT_ROOT="/app/Shared/Input"
OUTPUT_ROOT="/app/Shared/Output"

# Fixed file names inside each subfolder
CLINICAL_FILE_NAME="clinical.csv"
DEPMAP_FILE_NAME="depmap.csv"

SUBDIR_NAME="data_$(date +%Y-%m-%d_%H-%M-%S)"

INPUT_DIR="$INPUT_ROOT/$SUBDIR_NAME"
CLINICAL_FILE="$INPUT_DIR/$CLINICAL_FILE_NAME"
DEPMAP_FILE="$INPUT_DIR/$DEPMAP_FILE_NAME"

mkdir "$INPUT_DIR"

cp /app/Data/clin_data.csv "$CLINICAL_FILE"
cp "/app/Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv" "$DEPMAP_FILE"

# Infinite loop
while true; do

  for dir in "$OUTPUT_ROOT"/*; do
    # Skip if not a directory
    [ -d "$dir" ] || continue

    sleep 5 # Wait for output being ready

	rm -f "$dir"/*
	rmdir "$dir"

	SUBDIR_NAME="data_$(date +%Y-%m-%d_%H-%M-%S)"

    INPUT_DIR="$INPUT_ROOT/$SUBDIR_NAME"
    CLINICAL_FILE="$INPUT_DIR/$CLINICAL_FILE_NAME"
    DEPMAP_FILE="$INPUT_DIR/$DEPMAP_FILE_NAME"
	
	mkdir "$INPUT_DIR"
	
	cp /app/Data/clin_data.csv "$CLINICAL_FILE"
	cp "/app/Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv" "$DEPMAP_FILE"

  done

  # Wait before next scan
  sleep 5
done


