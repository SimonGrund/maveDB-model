library(xgboost)
library(caret)
library(tidyverse)
library(data.table)
library(ggpubr)

df = fread("Data/formatted_BRCA1_depletion_scores.tsv")
df$Score = df$Freq

run_boosted_model_comparison <- function(df, n_random = 1e3, seed = 42) {
  set.seed(seed)
  
  # Ensure factors
  df$WT <- as.factor(df$WT)
  df$MT <- as.factor(df$MT)
  
  # Prepare design matrix
  X <- model.matrix(Score ~ WT + MT + Location, data = df)[, -1]
  y_true <- df$Score
  
  # Train on true data
  dtrain_true <- xgb.DMatrix(data = X, label = y_true)
  
  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    max_depth = 8,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  
  model_true <- xgb.train(
    params = params,
    data = dtrain_true,
    nrounds = 100,
    verbose = 0
  )
  
  preds_true <- predict(model_true, dtrain_true)
  preds_true = pmin(pmax(preds_true, 0), 1) 
  rmse_true <- sqrt(mean((preds_true - y_true)^2))
  
  # Train on randomized Scores n times
  rmse_random <- numeric(n_random)
  
  for (i in 1:n_random) {
    y_shuffled <- sample(y_true)  # Shuffle Scores
    dtrain_rand <- xgb.DMatrix(data = X, label = y_shuffled)
    
    model_rand <- xgb.train(
      params = params,
      data = dtrain_rand,
      nrounds = 100,
      verbose = 0
    )
    
    preds_rand <- predict(model_rand, dtrain_rand)
    preds_rand = pmin(pmax(preds_rand, 0), 1) 
    
    rmse_random[i] <- sqrt(mean((preds_rand - y_shuffled)^2))
  }
  
  return(list(
    model_true = model_true,
    preds_true = preds_true,
    rmse_true = rmse_true,
    rmse_random = rmse_random,
    improvement = mean(rmse_random) - rmse_true
  ))
}

tt = run_boosted_model_comparison(df = df ) #Run the model

cat("RMSE on true data:", tt$rmse_true, "\n")
cat("Mean RMSE on randomized data:", mean(tt$rmse_random), "\n")
cat("Improvement over random models:", tt$improvement, "\n")

# plot feature importance
importance <- xgb.importance(model = tt$model_true)
importance_df <- as.data.frame(importance)
ggplot(importance_df, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Feature Importance", x = "Features", y = "Gain") +
  ggpubr::theme_classic2()

ggsave("Results/feature_importance_plot.pdf", width = 8, height = 14)

# Plot predictions vs true values
ggplot(data = data.frame(True = as.factor(round(df$Score,2)), Predicted = tt$preds_true), aes(x = True, y = Predicted)) +
 # geom_point(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA) +
#  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Predictions vs True Values", x = "True Values", y = "Predicted Values") +
  ggpubr::theme_classic2()

ggsave("Results/predictions_vs_true_values.pdf", width = 8, height = 6)

# Plot RMSE distribution
ggplot(data = data.frame(RMSE = tt$rmse_random), aes(x = RMSE)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  geom_vline(xintercept = tt$rmse_true, color = "red", linetype = "dashed") +
  labs(title = "Distribution of RMSE on Randomized Data", x = "RMSE", y = "Frequency") +
  ggpubr::theme_classic2()

ggsave("Results/rmse_distribution_plot.png", width = 8, height = 6)

# Save the model
xgb.save(tt$model_true, "Results/boosted_model.xgb")

#Save the list tt
saveRDS(tt, "Results/boosted_model_results.rds")
