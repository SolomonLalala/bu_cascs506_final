# import library
library(tidyverse)
library(knitr)
library(caret)
library(glmnet)
library(ggplot2)
library(GEOquery)
library(foreach)
library(doParallel)

#' Parse GEO dataset into expression matrix and phenotype table
#' @param dest_dir Destination directory to save parsed files
#' @param gse_id GEO Series ID for the dataset
#' @param transpose Logical, whether to transpose expression matrix
#' @param save Logical, whether to save output files
#' @return A list of expression dataframe and phenotype dataframe
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
    
    # save the expression df and phenotype df
    saveRDS(expr_df, file.path(dest_dir, paste0(gse_id, "_expr_df.rds")))
    saveRDS(pheno_df, file.path(dest_dir, paste0(gse_id, "_pheno_df.rds")))
  }
  
  # return the expression df and phenotype data frame
  return(list(expr_df, pheno_df))
}

#' Extract metadata and merge it with expression matrix
#'
#' @param gse_id GEO accession ID (e.g., "GSE40279")
#' @param dest_dir Directory where parsed RDS files are stored and merged result will be saved
#' @param save Logical; whether to save the merged dataset as RDS
#' @return A dataframe combining methylation expression and sample metadata

extract_merge <- function(gse_id, dest_dir, save = TRUE) {
  # load data
  expr <- readRDS(file.path(dest_dir, paste0(gse_id, "_expr_df.rds")))
  meta <- readRDS(file.path(dest_dir, paste0(gse_id, "_pheno_df.rds")))
  
  # extract relevant feature in metadata
  meta <- meta %>%
    select(
      sample_id = geo_accession,
      age = contains("age"),
      ethnicity = contains("ethnicity"),
      sex = contains("gender")
    ) %>%
    mutate(age = as.numeric(age))
  
  # add sample id to the expression dataframe
  expr$sample_id <- rownames(expr)
  
  # merge expression dataframe with meta
  merged <- inner_join(expr, meta, by = "sample_id")
  
  if (save) {
    # create output directory if it doesn't exist
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    
    saveRDS(merged, file = file.path(dest_dir, paste0(gse_id, "merged.rds")))
  }
  
  return(merged)
}


#' View a subset of the merged dataset
#' @param path Path to the RDS file that contains the merged dataset
#' @return A df with the first 6 samples and up to 10 columns
view_merged <- function(path) {
  
  # read in the data
  merged <- readRDS(path)
  
  merged_subset <- merged %>%
    select(
      sample_id,
      age,
      ethnicity,
      sex,
      everything()  # Then all other columns (CpGs)
    ) %>%
    slice(1:6) %>%
    select(1:10)
  
  return(merged_subset)
}


#' Stratified train-test split
#' @param dataset Input dataframe
#' @param train_frac Fraction of data to include in the training set
#' @param stratifier string of column name for stratification (e.g., age)
#' @param dest_dir Directory to save output files
#' @param save Logical, whether to save split datasets
#' @return A list containing train and test datasets
split_train_test <- function(
    dataset, train_frac, stratifier, dest_dir, save =TRUE
) 
{
  set.seed(123)
  
  dataset <- dataset %>%
    mutate(bin = ntile(.data[[stratifier]], 4))
  
  train_data <- dataset %>%
    group_by(bin, sex) %>%
    sample_frac(train_frac) %>%
    ungroup() %>%
    select(-bin)
  
  test_data <- anti_join(dataset, train_data, by = "sample_id") %>%
    select(-bin)
  
  data_list <- list(train_data, test_data)
  
  if (save) {
    # create output directory if it doesn't exist
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    
    
    saveRDS(train_data, file.path(dest_dir, "train_set.rds"))
    saveRDS(test_data, file.path(dest_dir, "test_set.rds"))
  }
  
  return(data_list)
  
}

#' Select features based on bootstrap frequency threshold
#' @param bootstrap_results Dataframe of CpGs selected in bootstraps
#' @param threshold Frequency ratio cutoff (e.g., 0.5)
#' @param dest_dir Directory to save selected features
#' @param save Logical, whether to save output
#' @return A vector of selected CpGs
select_feature <- function(bootstrap_results, threshold, dest_dir, save = TRUE) {
  feature_summary <- bootstrap_results %>%
    group_by(CpG) %>%
    summarize(Freq = n()) %>% # how many times each feature appears in bootstraps
    mutate(FreqRatio = Freq / length(unique(bootstrap_results$Bootstrap))) %>% # how often each feature appears in bootstraps
    arrange(desc(FreqRatio))
  
  selected_feature <- feature_summary %>%
    filter(FreqRatio > threshold) %>%
    pull(CpG)
  
  if (save) {
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    
    data_name <- deparse(substitute(bootstrap_results))
    filename <- paste0(data_name, "_selected_feature.rds")
    
    saveRDS(
      selected_feature,
      file.path(dest_dir, filename)
    )
  }
  
  return(selected_feature)
}

#' Evaluate model performance on test data
#' @param model Trained glmnet model
#' @param X_test Test set dataframe
#' @param y_test Test set labels
#' @param dest_dir Directory to save results
#' @param save Logical, whether to save evaluation metrics
#' @return A list containing predicted results and evaluation metrics
evaluate_model <- function(model, X_test, y_test, dest_dir, save = TRUE) {
  
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
  
  eval_list <- list(
    results_df = results_df,
    metrics = metrics
  )
  
  # save the results and metrics
  if (save) {
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    
    data_name <- deparse(substitute(model))
    filename <- paste0(data_name, "_res_and_eval.rds")
    
    saveRDS(
      eval_list,
      file.path(dest_dir, filename)
    )
  }
  
  return(eval_list)
}



#' Plot age distribution histogram
#' @param data Dataframe with an "age" column
#' @param dest_dir Directory to save the plot
#' @param save Logical, whether to save the plot
#' @return ggplot object
plot_age_distribution <- function(
    data, dest_dir, save = TRUE){
  
  # plot
  p_age <- ggplot(data = data, mapping = aes(x = age)) +
    geom_histogram(binwidth = 5, fill = "black", color = "white") +
    labs(title = "Distribution of Age",
         x = "Age (years)",
         y = "Count (# of sample)") +
    scale_x_continuous(breaks = seq(20, 100, by = 20)) +
    scale_y_continuous(breaks = seq(0, 100, by = 20)) +
    theme_minimal()
  
  # save the plot
  if (save) {
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    # --- Capture the name of the data input
    data_name <- deparse(substitute(data))
    filename <- paste0(data_name, "_age_distribution_plot.png")
    ggsave(
      filename = filename,
      plot = p_age, 
      path = dest_dir
    )
  }
  
  return(p_age)
}


#' Plot predicted vs actual age scatter plot
#' @param results Dataframe with Actual_Age and Predicted_Age
#' @param metrics List of evaluation metrics (correlation, RMSE et al)
#' @param dest_dir Directory to save the plot
#' @param save Logical, whether to save the plot
#' @return ggplot object
plot_evaluation <- function(
    results, metrics, dest_dir, save = TRUE
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
  
  
  # save the plot
  if (save) {
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE)
    }
    # --- Capture the name of the data input
    data_name <- deparse(substitute(results))
    filename <- paste0(data_name, "_age_prediction_plot.png")
    ggsave(
      filename = filename,
      plot = p, 
      path = dest_dir
    )}
  
  return(p)
  
}




