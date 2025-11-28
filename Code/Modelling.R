library(xgboost)
library(caret)
library(tidyverse)
library(data.table)
library(ggpubr)
library(doParallel)

set.seed(42)

# Enable parallel backend across all logical cores (hyperthreading)
available_cores <- parallel::detectCores(logical = TRUE)
if (!is.na(available_cores) && available_cores > 1) {
  cl <- parallel::makeCluster(available_cores)
  doParallel::registerDoParallel(cl)
  message(sprintf("Parallel backend registered with %d cores", available_cores))
} else {
  message("Parallel not enabled: single-core or detection failed")
}

# Load latest data
df <- fread("tmp/LATEST_formatted_data.tsv")

# Identify response column (prefer lower-case 'score')
response_col <- if ("score" %in% names(df)) {
  "score"
} else if ("Score" %in% names(df)) {
  "Score"
} else if ("RNAi" %in% names(df)) {
  # Backward compatibility if older file is present
  names(df)[names(df) == "RNAi"] <- "score"
  "score"
} else {
  stop("No response column named 'score'/'Score'/'RNAi' found in input data.")
}

# Drop columns with too many unique levels (likely IDs) and keep predictors
predictor_cols <- setdiff(names(df), response_col)

# Convert character columns to factor; drop very high-cardinality factors
sel_cols <- c(response_col, predictor_cols)
df_proc <- df[, ..sel_cols]
for (cn in predictor_cols) {
  if (is.character(df_proc[[cn]])) {
    df_proc[[cn]] <- as.factor(df_proc[[cn]])
  }
  if (is.factor(df_proc[[cn]])) {
    nlev <- nlevels(df_proc[[cn]])
    if (nlev > 0.5 * nrow(df_proc)) {
      df_proc[[cn]] <- NULL
      predictor_cols <- setdiff(predictor_cols, cn)
    }
  }
}

# Split train/test
train_idx <- caret::createDataPartition(df_proc[[response_col]], p = 0.8, list = FALSE)
train_df <- df_proc[train_idx]
test_df  <- df_proc[-train_idx]

# Use caret dummyVars to one-hot encode predictors consistently
predictor_names <- setdiff(names(train_df), response_col)
X_df_train <- as.data.frame(train_df[, ..predictor_names])
X_df_test  <- as.data.frame(test_df[, ..predictor_names])

dummies <- caret::dummyVars(~ ., data = X_df_train, fullRank = TRUE)
X_train <- predict(dummies, newdata = X_df_train)
X_test  <- predict(dummies, newdata = X_df_test)

y_train <- train_df[[response_col]]
y_test  <- test_df[[response_col]]

# Ensure numeric response
y_train <- as.numeric(y_train)
y_test  <- as.numeric(y_test)

# Scale target to [-1, 1] using train min/max
scale_to_pm1 <- function(vec, mn, mx) {
  if (is.na(mn) || is.na(mx) || mx == mn) {
    return(rep(0, length(vec)))
  }
  -1 + 2 * (vec - mn) / (mx - mn)
}
train_min <- min(y_train, na.rm = TRUE)
train_max <- max(y_train, na.rm = TRUE)
y_train <- scale_to_pm1(y_train, train_min, train_max)
y_test  <- scale_to_pm1(y_test,  train_min, train_max)

# Drop rows with NA targets in train
keep_idx <- stats::complete.cases(y_train)
if (!all(keep_idx)) {
  X_train <- X_train[keep_idx, , drop = FALSE]
  y_train <- y_train[keep_idx]
}

cat(sprintf("Dims before NZV: X_train=%d x %d, y_train=%d, X_test=%d x %d\n",
            nrow(X_train), ncol(X_train), length(y_train), nrow(X_test), ncol(X_test)))

# Remove near-zero variance predictors to speed up training
nzv <- caret::nearZeroVar(X_train)
if (length(nzv) > 0) {
  X_train <- X_train[, -nzv, drop = FALSE]
  X_test  <- X_test[, -nzv, drop = FALSE]
}

cat(sprintf("Dims after NZV: X_train=%d x %d, y_train=%d, X_test=%d x %d\n",
            nrow(X_train), ncol(X_train), length(y_train), nrow(X_test), ncol(X_test)))

# Caret training with a compact grid and repeated CV
tr_ctrl <- caret::trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 1,
  verboseIter = FALSE,
  allowParallel = TRUE
)

grid <- expand.grid(
  nrounds = c(300, 600),            # allow more trees for lower eta
  max_depth = c(3, 4, 6),           # shallower trees for regularization
  eta = c(0.01, 0.03, 0.05, 0.1),   # stronger shrinkage options
  gamma = c(0, 0.5, 1),             # min loss reduction for split
  colsample_bytree = c(0.6, 0.8),   # feature subsampling
  min_child_weight = c(3, 5, 7),    # discourage overly specific splits
  subsample = c(0.6, 0.8)           # row subsampling
)

message("Training XGBoost model with cross-validation...")
fit <- caret::train(
  x = as.data.frame(X_train),
  y = y_train,
  method = "xgbTree",
  trControl = tr_ctrl,
  tuneGrid = grid,
  metric = "RMSE",
  nthread = available_cores
)

# Evaluate on test set
pred_test <- predict(fit, newdata = as.data.frame(X_test))
rmse_test <- sqrt(mean((pred_test - y_test)^2))
mae_test  <- mean(abs(pred_test - y_test))
r2_test   <- 1 - sum((pred_test - y_test)^2) / sum((y_test - mean(y_test))^2)
spearman  <- suppressWarnings(cor(pred_test, y_test, method = "spearman", use = "complete.obs"))

cat(sprintf("Test RMSE: %.4f\n", rmse_test))
cat(sprintf("Test MAE: %.4f\n", mae_test))
cat(sprintf("Test R2: %.4f\n", r2_test))
cat(sprintf("Test Spearman: %.4f\n", spearman))

# Shuffle-baseline: train with shuffled y_train, evaluate on true y_test
n_random <- 100
best <- fit$bestTune

params_best <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  max_depth = best$max_depth,
  eta = best$eta,
  gamma = best$gamma,
  colsample_bytree = best$colsample_bytree,
  min_child_weight = best$min_child_weight,
  subsample = best$subsample,
  nthread = available_cores
)

dtest  <- xgb.DMatrix(data = X_test)
rmse_random <- numeric(n_random)

message("Running shuffled-label baseline...")
for (i in seq_len(n_random)) {
  y_shuf <- sample(y_train)
  dtrain_shuf <- xgb.DMatrix(data = X_train, label = y_shuf)
  mod_shuf <- xgb.train(params = params_best, data = dtrain_shuf, nrounds = best$nrounds, verbose = 0)
  pred_shuf <- predict(mod_shuf, dtest)
  rmse_random[i] <- sqrt(mean((pred_shuf - y_test)^2))
}

improvement <- mean(rmse_random) - rmse_test
cat(sprintf("Mean shuffled RMSE: %.4f\n", mean(rmse_random)))
cat(sprintf("Improvement over shuffled: %.4f\n", improvement))

# Feature importance (gain)
final_booster <- fit$finalModel
imp <- xgb.importance(model = final_booster, feature_names = colnames(X_train))
imp_df <- as.data.frame(imp)

dir.create("Results", showWarnings = FALSE)
fwrite(imp_df, file = "Results/feature_importance_gain.csv")

# Plot top 30 important features
top_n <- 30
var_train <- apply(X_train, 2, var, na.rm = TRUE)
threshold_var <- quantile(var_train, 0.30)
keep_var <- which(var_train > threshold_var)
X_train <- X_train[, keep_var, drop = FALSE]
X_test  <- X_test[, keep_var, drop = FALSE]
cat(sprintf("Dims after variance filter (30%% low removed): X_train=%d x %d, X_test=%d x %d\n",
            nrow(X_train), ncol(X_train), nrow(X_test), ncol(X_test)))
imp_plot_df <- imp_df %>% dplyr::slice(seq_len(min(nrow(imp_df), top_n))) %>%
  dplyr::mutate(Feature = factor(Feature, levels = rev(Feature)))

ggplot(imp_plot_df, aes(x = Feature, y = Gain)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Feature Importances (Gain)", x = "Feature", y = "Gain") +
  ggpubr::theme_classic2()
ggsave("Results/feature_importance_top30.pdf", width = 8, height = 10)

# Predictions vs True (test set)
ggplot(data.frame(True = y_test, Predicted = pred_test), aes(x = True, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Predicted vs True (Test)", x = "True score", y = "Predicted score") +
  ggpubr::theme_classic2()
ggsave("Results/pred_vs_true_test.pdf", width = 7, height = 6)

# RMSE distribution for shuffled baseline
ggplot(data.frame(RMSE = rmse_random), aes(x = RMSE)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = rmse_test, color = "red", linetype = "dashed") +
  labs(title = "RMSE on Shuffled-Label Models (Test)", x = "RMSE", y = "Count") +
  ggpubr::theme_classic2()
ggsave("Results/rmse_random_distribution.pdf", width = 7, height = 6)

# Save model artifacts and metrics
saveRDS(fit, file = "Results/xgb_caret_model.rds")
try({ xgb.save(final_booster, fname = "Results/xgb_booster.model") }, silent = TRUE)

metrics <- data.frame(
  rmse_test = rmse_test,
  mae_test = mae_test,
  r2_test = r2_test,
  spearman_test = spearman,
  rmse_random_mean = mean(rmse_random),
  rmse_random_sd = sd(rmse_random),
  improvement = improvement
)
fwrite(metrics, file = "Results/metrics_test.csv")
fwrite(data.table(rmse_random = rmse_random), file = "Results/rmse_random_values.csv")

message("Done. Artifacts saved in Results/")

# Clean up parallel backend if started
try({
  if (exists("cl") && inherits(cl, "cluster")) parallel::stopCluster(cl)
}, silent = TRUE)
