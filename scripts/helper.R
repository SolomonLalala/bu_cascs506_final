# import library
library(tidyverse)
library(caret)
library(glmnet)
library(ggplot2)
library(GEOquery)
library(foreach)
library(doParallel)

# query and parse data from the GEO database
parse_geo <- function(dest_dir, gse_id, transpose = TRUE, save = TRUE) {
  
  # query geo and get the gse object
  gse_obj <- getGEO(
    GEO = gse_id, 
    destdir = dest_dir, # Specify the destination directory
    GSEMatrix = TRUE, # download the GSE Series Matrix
    getGPL = FALSE, # do not download the GPL file
    parseCharacteristics = TRUE # parse characteristics, default is TRUE
  )
  
  # get expression set object, expression matrix, and phenotype dataframe
  eset <- gse_obj[[1]]
  expr_mtx <- exprs(eset)
  pheno_df <- pData(eset)
  
  # transpose the expression matrix
  if (transpose) {
    expr_mtx_t <- t(expr_mtx)
    expr_df <- as.data.frame(expr_mtx_t)
  } else {
    expr_df <- as.data.frame(expr_mtx)
  }
  
  if (save) {
    # create output directory if it doesn't exist
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    
    # save the expression df and phenotype data frame to CSV files
    saveRDS(expr_df, file.path(dest_dir, paste0(gse_id, "_expr_df.rds")))
    saveRDS(pheno_df, file.path(dest_dir, paste0(gse_id, "_pheno_df.rds")))
  }
  
  # return the expression df and phenotype data frame
  return(list(expr_df, pheno_df))
}


split_train_test <- function(
    dataset, train_frac, stratifier, dest_dir, save =TRUE
) 
{
  set.seed(123)
  
  stratifier1 = stratifier[[1]]
  
  dataset <- dataset %>%
    mutate(bin = ntile(.data[[stratifier1]], 4))
  
  train_data <- dataset %>%
    group_by(bin, sex) %>%
    sample_frac(train_frac) %>%
    ungroup() %>%
    select(-bin)
  
  test_data <- anti_join(dataset, train_data, by = "sample_id") %>%
    select(-bin)
  
  if (save) {
    # create output directory if it doesn't exist
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    
    saveRDS(train_data, file.path(dest_dir, "train_set.rds"))
    saveRDS(test_data, file.path(dest_dir, "test_set.rds"))
  }
  
  return(list(train_data, test_data))
  
}

select_feature <- function(bootstrap_results, threshold, dest_dir) {
  feature_summary <- bootstrap_results %>%
    group_by(CpG) %>%
    summarize(Freq = n()) %>% # how many times each feature appears in bootstraps
    mutate(FreqRatio = Freq / length(unique(bootstrap_results$Bootstrap))) %>% # how often each feature appears in bootstraps
    arrange(desc(FreqRatio))
  
  selected_feature <- feature_summary %>%
    filter(FreqRatio > threshold) %>%
    pull(CpG)
  
  saveRDS(
    selected_feature,
    file.path(dest_dir, "selected_feature.rds")
  )
  
  return(selected_feature)
}

evaluate_model <- function(model, X_test, y_test, dest_dir) {
  
  y_pred <- predict(model, newx = as.matrix(X_test), s = lambda_min)
  y_pred <- as.numeric(y_pred)
  results_df <- data.frame(
    Sample = test_set$sample_id,
    Actual_Age = y_test,
    Predicted_Age = y_pred
  )
  
  # Residual Sum of Squares: how far off the predicted values are from the actual values
  rss <- sum((y_test - y_pred)^2) 
  # Total Sum of Squares: how far off the actual values are from their mean
  tss <- sum((y_test - mean(y_test))^2) 
  # coefficient of determination: how much of the variance in the data is explained by the model
  r_squared <- 1 - rss/tss 
  # correlation coeffcient: how well the predicted values correlate with the actual values
  r <- cor(y_test, y_pred)
  # Root Mean Squared Error: average of the square root of the difference between the predicted and actual values
  rmse <- sqrt(mean((y_test - y_pred)^2))
  
  # wrap the metrics in a list
  metrics <- list(
    rss = rss,
    tss = tss,
    r_squared = r_squared,
    r = r,
    rmse = rmse
  )
  
  # save the results and metrics
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  saveRDS(
    results_df,
    file.path(dest_dir, "results_df.rds")
  )
  saveRDS(
    metrics,
    file.path(dest_dir, "evaluation_metrics.rds")
  )
  
  eval_list <- list(
    results_df = results_df,
    metrics = metrics
  )
  
  # return 
  return(eval_list)
}

plot_evaluation <- function(
    results, metrics, dest_dir
){
  
  r <- metrics$r
  rmse <- metrics$rmse
  
  # fit_line <- lm(y_pred ~ y_test)
  # slope_fit <- coef(fit_line)[2]
  # intercept_fit <- coef(fit_line)[1]
  
  p <- ggplot(results, aes(x = Actual_Age, y = Predicted_Age)) +
    geom_point() +
    geom_abline(
      slope = 1, 
      intercept = 0, 
      linetype = "dashed", 
      color = "red") +
    # geom_abline(
    #   slope = slope_fit, 
    #   intercept = intercept_fit, 
    #   linetype = "solid",
    #   color = "blue") +
    labs(title = "Predicted vs Actual Age",
         x = "Actual Age (years)",
         y = "Predicted Age (years)") +
    annotate("text", x = min(y_test), y = max(y_test), 
             label = paste0("r = ", round(r, 2), 
                            "\nrmse = ", round(rmse, 2)),
             hjust = 0, vjust = 1, size = 4) +
    scale_x_continuous(limits = c(0, 100)) + 
    scale_y_continuous(limits = c(0, 100)) + 
    theme_minimal()
  
  ggsave(
    filename = "age_prediction_plot.png",
    plot = p, 
    path = dest_dir
    )
  
  return(p)
}




