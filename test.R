# setup
library(testthat)
source("scripts/helper.R")
data_dir <- "data" 
fig_dir <- "figures_test" 
gse_accession <- "GSE40279" 

test_that("Required libraries are loaded", {
  required_libs <- c("tidyverse", "knitr", "caret", "glmnet", "ggplot2", "GEOquery", "foreach", "doParallel", "testthat")
  for (pkg in required_libs) {
    expect_true(pkg %in% loadedNamespaces(), info = paste(pkg, "is not loaded"))
  }
})

# load data
merged <- readRDS(file.path(data_dir, paste0(gse_accession, "merged.rds")))
test_that("Merged dataset loads with expected dimensions", {
  expect_true(ncol(merged) > 10)
  expect_true(nrow(merged) >= 6)
})

p_age_merged <- plot_age_distribution(
  data = merged,
  dest_dir = fig_dir, 
  save = TRUE)
p_age_merged

# split data into train and test sets
train_test_list <- split_train_test(
  dataset = merged, 
  train_frac = 0.75, 
  stratifier = "age", 
  dest_dir = data_dir, 
  save = TRUE)
train_set <- train_test_list[[1]]
test_set <- train_test_list[[2]]

X_train <- train_set %>%
  select(-sample_id, -age, -ethnicity, -sex) 
y_train <- train_set$age
X_test <- test_set %>%
  select(-sample_id, -age, -ethnicity, -sex)
y_test <- test_set$age
saveRDS(X_train, file = file.path(data_dir, "X_train.rds"))
saveRDS(y_train, file = file.path(data_dir, "y_train.rds"))
saveRDS(X_test, file = file.path(data_dir, "X_test.rds"))
saveRDS(y_test, file = file.path(data_dir, "y_test.rds"))

# Standardize
# obtain mean and sd
preProcValues <- preProcess(as.matrix(X_train), method = c("center", "scale"))
# Apply the transformation to the data
X_train_standardized <- predict(preProcValues, as.matrix(X_train))
X_test_standardized <- predict(preProcValues, as.matrix(X_test))
saveRDS(preProcValues, file = file.path(data_dir, "preProcValues.rds"))
saveRDS(X_train_standardized, file = file.path(data_dir, "X_train_standardized.rds"))
saveRDS(X_test_standardized, file = file.path(data_dir, "X_test_standardized.rds"))

test_that("Standardized data has correct shape", {
  expect_equal(dim(X_train_standardized)[2], ncol(X_train))
  expect_equal(dim(X_test_standardized)[2], ncol(X_test))
})

# Train model
alpha = 0.5
n_cores <- parallel::detectCores() - 1
n_fold = 5
iterations = 100
model_params <- list(
  n_cores = n_cores,
  alpha = alpha,
  n_fold = n_fold,
  iterations = iterations
)
saveRDS(model_params, file = file.path(data_dir, "model_params.rds"))

# cv
source("scripts/cv.R")
cv_model <- readRDS(file.path(data_dir, "cv_model.rds"))
plot(cv_model)
lambda_min <- cv_model$lambda.min
lambda_min
ggsave(
  filename = file.path(fig_dir, "cv_model_lambda.png"),
  plot = plot(cv_model),
)
saveRDS(lambda_min, file = file.path(data_dir, "lambda_min.rds"))

# Test cv.glmnet result
test_that("cv.glmnet returns valid cv model with lambda.min", {
  expect_s3_class(cv_model, "cv.glmnet")
  expect_true(!is.null(cv_model$lambda.min))
  expect_true(is.numeric(cv_model$lambda.min))
  expect_gt(cv_model$lambda.min, 0)
})

# bootstrap
source("scripts/bootstrap.R")

bootstrap_res <- readRDS(file.path(data_dir, "bootstrap_results.rds"))
# Test bootstrap result structure
test_that("run_bootstrap returns non-empty bootstrap matrix", {
  expect_true(is.data.frame(bootstrap_res))
  expect_gt(ncol(bootstrap_res), 1)
  expect_gt(nrow(bootstrap_res), 1)
  expect_true(all(bootstrap_res$Bootstrap >= 1 & bootstrap_res$Bootstrap <= iterations))  # frequency-like values
})

# Feature selection
selected_feature <- select_feature(
  bootstrap_results = bootstrap_res, 
  threshold = 0.1,
  dest_dir = data_dir,
  save = TRUE)
# Test select_feature result
test_that("select_feature returns feature names", {
  expect_true(is.character(selected_feature))
  expect_gte(length(selected_feature), 1)
  expect_true(all(selected_feature %in% colnames(X_train)))
})


# Train final model
X_train_standardized_final <- as.matrix(as.data.frame(X_train_standardized)[, selected_feature])
X_test_standardized_final <- as.matrix(as.data.frame(X_test_standardized)[, selected_feature])

set.seed(123)
final_model <- glmnet(
  X_train_standardized_final, 
  y_train, 
  alpha = alpha, 
  lambda = lambda_min, 
  type.measure = "mse", 
  family = "gaussian", 
  standardize = FALSE) 

saveRDS(final_model, file = file.path(data_dir, "final_model.rds"))

# Test glmnet final model
test_that("final_model is a valid", {
  expect_s3_class(final_model, "glmnet")})

# evaluation
evaluation <- evaluate_model(
  model = final_model, 
  X_test = X_test_standardized_final, 
  y_test = y_test, 
  dest_dir = data_dir)

eval_res <- evaluation$results_df
eval_metrics <- evaluation$metrics

test_that("Evaluation returns expected metrics", {
  expect_equal(names(eval_metrics)[4], "r")  # check the 4th name
  
  # Check that the value of eval_metrics[[4]] is numeric
  expect_true(is.numeric(eval_metrics[[4]]))
})

plot_eval <- plot_evaluation(
  results = eval_res, 
  metrics = eval_metrics, 
  dest_dir = fig_dir)
plot_eval

sessionInfo()





