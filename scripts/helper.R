# import library
library(tidyverse)
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
    
    write_rds(train_data, file.path(dest_dir, "train_data.rds"))
    write_rds(test_data, file.path(dest_dir, "test_data.rds"))
  }
  
  return(list(train_data, test_data))
  
}


select_features <- function(bootstrap_results, threshold) {
  cpg_summary <- bootstrap_results %>%
    group_by(CpG) %>%
    summarize(Freq = n()) %>%
    mutate(FreqRatio = Freq / length(unique(bootstrap_results$Bootstrap))) %>%
    arrange(desc(FreqRatio))
  
  selected_cpgs <- cpg_summary %>%
    filter(FreqRatio > 0.5) %>%
    pull(CpG)
  
  saveRDS(
    selected_cpgs, 
    file = "data/selected_cpgs.rds"
  )
  
  return(selected_cpgs)
}

train_model <- function(X_train, y_train, alpha, lambda, features) {
  X_train <- as.matrix(X_train[, features])
  y_train <- X_train$age
  
  X_train_scaled <- scale(X_train)
  
  best_lambda <- readRDS("data/cv_model.rds")$lambda.min
  
  model <- glmnet(
    X_train_scaled, 
    y_train, 
    alpha = 0.5, 
    lambda = best_lambda, 
    standardize = FALSE)
}



