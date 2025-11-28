# Pipeline runner for maveDB-model
# Description of clinical input format:
# - Two columns CSV/TSV
# - Column 1: depmap_id (string identifier matching DepMap IDs)
# - Column 2: RNAi (numeric score, outcome to predict)
#
# Configure input paths here or pass via command line args.
# Example CLI:
#   Rscript Code/run_pipeline.R Data/clin_data.csv Data/depmap_export_2025-11-26\ 09_50_16.184205_subsetted.csv

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
clin_path <- if (length(args) >= 1) args[1] else "Data/clin_data.csv"
depmap_path <- if (length(args) >= 2) args[2] else "Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv"

message(sprintf("Using clin_data: %s", clin_path))
message(sprintf("Using depmap export: %s", depmap_path))

# Expose paths to data preparation script via global variables
assign("CLIN_DATA_PATH", clin_path, envir = .GlobalEnv)
assign("DEPMAP_EXPORT_PATH", depmap_path, envir = .GlobalEnv)

# 1) Load and format data
message("[1/3] Running data load and formatting...")
source("Code/Load and format data.R")

# 2) Train model and produce artifacts
message("[2/3] Training model and generating artifacts...")
source("Code/Modelling.R")

# 3) Render report (HTML only) and copy HTML to Results
message("[3/3] Rendering report (HTML)...")
qmd <- "Report/model_report.qmd"
if (!requireNamespace("quarto", quietly = TRUE)) {
  message("Quarto not available as an R package; invoking CLI if present...")
}
# Prefer CLI for broader compatibility
qmd_cmd_html <- sprintf("quarto render %s --to html", shQuote(qmd))

# Try rendering via system() calls (HTML only)
html_status <- system(qmd_cmd_html)
if (html_status == 0 && file.exists("Report/model_report.html")) {
  dir.create("Results", showWarnings = FALSE)
  file.copy("Report/model_report.html", "Results/model_report.html", overwrite = TRUE)
  message("Report HTML copied to Results/model_report.html")
} else {
  message("HTML render failed or file missing; please ensure Quarto is installed.")
}

message("Pipeline complete.")
