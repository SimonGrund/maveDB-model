# maveDB-model

Predicting BRCA1 variant functional impact using MaveDB depletion scores and XGBoost.

## Overview

This project builds a gradient-boosted regression model (XGBoost) to predict the functional impact of missense variants in the BRCA1 RING and BRCT domains, using multiplexed assay data from [MaveDB](https://mavedb.org/score-sets/urn:mavedb:00000081-a-1). The goal is to learn sequence-level determinants of variant pathogenicity from experimentally measured depletion scores.

## Datasets

### MaveDB BRCA1 Depletion Scores

- **Source:** [urn:mavedb:00000081-a-1](https://mavedb.org/score-sets/urn:mavedb:00000081-a-1)
- **File:** `Data/urn_mavedb_00000081-a-1_scores.csv`
- **Description:** Depletion scores for missense variants in the BRCA1 RING and BRCT domains. The score represents the number of replicates in which a variant was depleted relative to a control siRNA replicate. Each variant is described in HGVS protein notation (e.g. `p.Thr167Cys`).

### DepMap Cell Line Data

- **Source:** [DepMap Portal](https://depmap.org/)
- **File:** `Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv`
- **Description:** Multi-omic features for breast cancer cell lines, including:
  - Gene expression (Expression Public 25Q3) for ~480 cancer-related genes
  - Absolute copy number (Omics Absolute CN Gene Public 24Q4) for the same gene set
  - Damaging mutation status (binary) across ~480 genes
  - Cell line metadata: lineage, subtype (e.g. HER2+), and DepMap ID

### PARP1 Dependency Scores

- **File:** `Data/clin_data.csv`
- **Description:** PARP1 gene effect (dependency) scores per cell line from DepMap CRISPR screens.

## Methodology

### 1. Data Preprocessing (`Code/Load and format data.R`)

- Parses HGVS protein notation into three features:
  - **WT** — wild-type amino acid (3-letter code)
  - **Location** — residue position in the protein
  - **MT** — mutant amino acid (3-letter code)
- Computes a normalised depletion frequency: `Freq = Score / Replicates`, capped at 1.
- Removes rows with missing values.

### 2. Modelling (`Code/Modelling.R`)

- **Algorithm:** XGBoost regression (`reg:squarederror`)
- **Features:** One-hot encoded WT amino acid, MT amino acid, and numeric residue position.
- **Target:** Normalised depletion frequency (0–1).
- **Hyperparameters:** max_depth = 8, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8, nrounds = 100.
- **Validation:** Permutation-based null comparison — the model is re-trained 1,000 times on shuffled labels, and the RMSE of the true model is compared against the null RMSE distribution.
- **Outputs:**
  - Trained model (`Results/boosted_model.xgb`)
  - Full results object (`Results/boosted_model_results.rds`)
  - Feature importance plot (`Results/feature_importance_plot.pdf`)
  - Predictions vs. true values (`Results/predictions_vs_true_values.pdf`)
  - Null RMSE distribution (`Results/rmse_distribution_plot.png`)

## About XGBoost

[XGBoost](https://xgboost.readthedocs.io/) (eXtreme Gradient Boosting) is an ensemble learning method that builds a sequence of decision trees, where each new tree corrects the errors of the previous ones. Key concepts:

- **Gradient boosting:** Trees are added iteratively, each one fitting the residual errors (negative gradient of the loss function) from the current ensemble. This additive strategy progressively reduces prediction error.
- **Regularisation:** XGBoost includes L1 and L2 penalties on leaf weights and limits tree complexity (`max_depth`), which reduces overfitting compared to traditional gradient boosting.
- **Stochastic sampling:** At each boosting round, only a random subset of rows (`subsample = 0.8`) and features (`colsample_bytree = 0.8`) are used, improving generalisation and reducing variance.
- **Learning rate (eta):** Controls the contribution of each tree. A smaller value (here `eta = 0.1`) requires more rounds but produces a smoother, more robust model.
- **Squared error objective:** The model minimises mean squared error between predicted and observed depletion frequencies, making it suitable for continuous targets bounded in [0, 1].

In this project, the input features (wild-type residue, mutant residue, and position) are one-hot encoded into a sparse design matrix. XGBoost natively handles sparse inputs efficiently, making it well-suited for categorical amino acid features expanded via `model.matrix()`.

## Project Structure

```
├── Code/
│   ├── Install_packages.R        # Dependency installation
│   ├── Load and format data.R    # Data parsing and feature engineering
│   └── Modelling.R               # XGBoost training and evaluation
├── Data/                          # Input datasets
├── Results/                       # Model outputs and plots
├── Docker/                        # Containerised environment
│   ├── Dockerfile
│   ├── build_image.sh
│   └── run-container.sh
├── run.sh                         # Entry point (runs Modelling.R)
└── README.md
```

## Requirements

R packages: `xgboost`, `caret`, `tidyverse`, `data.table`, `ggpubr`

Install via:

```r
source("Code/Install_packages.R")
```

## Running

```bash
# Locally
Rscript Code/Modelling.R

# Or via the wrapper script
bash run.sh

# Or via Docker
cd Docker && bash build_image.sh && bash run-container.sh
```
