library(glmnet)
library(foreach)
library(doParallel)
library(tidyverse)

train_model_bootstrap <- function(
    X_train, y_train, alpha, n_cores, iteration, lambda) {
  
  registerDoParallel(cores = n_cores)
  
  set.seed(123)
  
  bootstrap_results <- foreach(
    i = 1:iteration, .combine = rbind, .packages = "glmnet"
  ) %dopar% {
    boot_indices <- sample(1:nrow(X_train), replace = TRUE)
    X_boot <- X_train[boot_indices, ]
    y_boot <- y_train[boot_indices]
    
    model_boot <- glmnet(
      X_boot, y_boot, alpha = alpha, lambda = lambda, standardize = TRUE)
    coef_boot <- coef(model_boot)
    nonzero_cpgs <- rownames(coef_boot)[which(coef_boot != 0)]
    nonzero_cpgs <- nonzero_cpgs[nonzero_cpgs != "(Intercept)"]
    
    data.frame(CpG = nonzero_cpgs, Bootstrap = i)
  }
  
  # save results
  saveRDS(bootstrap_results, "data/bootstrap_results.rds")
  
  return(bootstrap_results)
}

# import data
X_train <- readRDS("data/X_train.rds")
y_train <- readRDS("data/y_train.rds")
model_params <- readRDS("data/model_params.rds")
# import lambda
lambda_min <- readRDS("data/lambda_min.rds")

X_train <- as.matrix(X_train)
alpha <- model_params$alpha
n_cores <- model_params$n_cores
n_fold <- model_params$n_fold
iterations <- model_params$iterations

bootstrap_res <- train_model_bootstrap(
  X_train = X_train, 
  y_train = y_train, 
  alpha = alpha, 
  n_cores = n_cores, 
  iteration = iterations, 
  lambda = lambda_min 
)

