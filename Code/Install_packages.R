# Define required packages
required_packages <- c(
  "data.table",
  "tidyverse",     # includes dplyr, tidyr, readr, ggplot2, etc.
  "xgboost",
  "RhpcBLASctl",
  "caret",
  "ggpubr",
  "doParallel",    # parallel backend for caret/xgboost
  "kableExtra"      # PDF table formatting in report
)

# Install any packages that are not already installed
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
  }
}