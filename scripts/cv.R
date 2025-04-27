library(glmnet)
library(foreach)
library(doParallel)
library(tidyverse)

train_model_cross_validation <- function(
    X_train, y_train, alpha, n_fold, n_cores) {
  
  registerDoParallel(cores = n_cores)
  
  set.seed(123)
  
  # cross validation
  cv_model <- cv.glmnet(
    x = X_train, 
    y = y_train,
    alpha = alpha,
    type.measure = "mse",
    nfolds = n_fold,
    family = "gaussian",
    standardize = TRUE
  )

  # save results
  saveRDS(cv_model, "data/cv_model.rds")

  return(cv_model)
}

# import data
X_train <- readRDS("data/X_train.rds")
y_train <- readRDS("data/y_train.rds")
model_params <- readRDS("data/model_params.rds")

X_train <- as.matrix(X_train)
alpha <- model_params$alpha
n_cores <- model_params$n_cores
n_fold <- model_params$n_fold

elastic_net_model <- train_model_cross_validation(
  X_train, 
  y_train, 
  alpha = alpha, 
  n_fold = n_fold, 
  n_cores = n_cores
  )
