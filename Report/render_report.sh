#!/usr/bin/env zsh
set -euo pipefail
cd "$(dirname $0)/.."  # go to project root
REPORT="Report/model_report.qmd"
HTML_OUT="Report/model_report.html"
PDF_OUT="Report/model_report.pdf"
RESULTS_PDF="Results/model_report.pdf"
RESULTS_HTML="Results/model_report.html"

if ! command -v quarto >/dev/null 2>&1; then
  echo "Quarto not found. Install with: brew install quarto" >&2
  exit 1
fi

echo "Rendering HTML..."
quarto render "$REPORT" --to html

echo "Rendering PDF..."
quarto render "$REPORT" --to pdf

if [ -f "$PDF_OUT" ]; then
  mkdir -p Results
  cp "$PDF_OUT" "$RESULTS_PDF"
  echo "Copied PDF to $RESULTS_PDF"
else
  echo "PDF not generated (missing $PDF_OUT)" >&2
  exit 2
fi

# Copy HTML report as well
if [ -f "$HTML_OUT" ]; then
  mkdir -p Results
  cp "$HTML_OUT" "$RESULTS_HTML"
  echo "Copied HTML to $RESULTS_HTML"
else
  echo "HTML not generated (missing $HTML_OUT)" >&2
fi

echo "Done. Outputs: $HTML_OUT and $RESULTS_PDF"