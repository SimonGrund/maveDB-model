# DepMap Modeling Pipeline (XGBoost + HTML Report)

End-to-end pipeline to train a predictive model on DepMap data and produce an HTML summary report. The pipeline ingests DepMap features (expression, copy-number, mutation status) and a simple clinical/outcome file, then trains an XGBoost model with cross-validation, evaluates against a shuffled-label baseline, exports feature importance, and renders a report.

## Inputs

- Clinical/outcome file: CSV/TSV with exactly two columns
  - Column 1: `depmap_id` (matches DepMap IDs)
  - Column 2: `RNAi` (numeric score, the outcome to predict)
- DepMap export CSV: feature matrix (expression/CNV/mutation columns). Example path:
  - `Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv`

## Outputs (in `Results/`)

- `metrics_test.csv`: RMSE, MAE, R2, Spearman, shuffled RMSE mean/sd, improvement
- `feature_importance_gain.csv`: XGBoost Gain importance for all features
- `rmse_random_values.csv`: RMSE distribution from shuffled-label models
- Figures: `pred_vs_true_test.pdf`, `rmse_random_distribution.pdf`, `feature_importance_top30.pdf`
- Report: `Results/model_report.html` (also `Report/model_report.html`)

Notes:

- Target is scaled to [-1, 1] based on the training set to avoid leakage.
- Parallelism: uses all logical cores via `doParallel` and sets xgboost `nthread`.

## Quick Start (Local)

Prerequisites: R (4.3+), Quarto (for report), and R packages.

1. Install R packages (first run only):

```zsh
cd /Users/simon/Documents/Projects/VUSdx/Mave-model/maveDB-model
Rscript Code/Install_packages.R
```

2. Run the full pipeline (defaults shown):

```zsh
./run.sh \
	Data/clin_data.csv \
	"Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv"
```

3. Render the report manually (HTML-only) if desired:

```zsh
cd Report
./render_report.sh
```

## Docker

This repository includes a Dockerfile and `docker-compose.yml` for reproducible runs.

- Build and run with compose:

```zsh
cd /Users/simon/Documents/Projects/VUSdx/Mave-model/maveDB-model
docker compose up --build
```

- Override input arguments (maps to `run.sh` args):

```zsh
docker compose run --rm mavedb-model \
	/app/Data/clin_data.csv \
	"/app/Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv"
```

Volumes mounted by compose:

- `./Results:/app/Results` (artifacts are persisted to the host)
- `./Data:/app/Data`
- `./tmp:/app/tmp`

By default the container entrypoint is `run.sh`, which invokes `Code/run_pipeline.R` with the provided args and copies the HTML report into `Results/`.

## Pipeline Details

- Data prep: `Code/Load and format data.R` loads and joins DepMap features with the provided clinical/outcome file and writes `tmp/LATEST_formatted_data.tsv`.
- Modeling: `Code/Modelling.R`
  - One-hot encoding via `caret::dummyVars`
  - Remove near-zero variance predictors; apply a variance filter
  - XGBoost via `caret::train(method = "xgbTree")` with CV grid
  - Shuffled-label baseline (n=100) to quantify improvement over chance
  - Feature importance (Gain) exported to CSV and plotted
- Report: `Report/model_report.qmd` (HTML only). The `render_report.sh` copies the HTML into `Results/`.

## Troubleshooting

- Report render requires Quarto on the host (`brew install quarto` on macOS). The Docker image can render HTML internally if Quarto is installed there; PDFs have been disabled by default.
- If Docker on Apple Silicon has base image issues, ensure you’re on a multi-arch base (we use `rocker/r-ver:4.3.2`) or specify a platform flag.

## License

This repository contains project-specific code. Unless specified otherwise, all rights reserved by the project owner.
