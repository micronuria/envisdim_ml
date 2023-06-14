
# Author: Nuria Fernández-Gonzalez

###########################################
#
# ASVs RF MODEL 1
#
###########################################

#---------------------------------------------
# Description. This script trains a RF model
# using 
# - ASVs and some environmental variables
# - Normalization method: proportions
# - Not filtered by prevalence
# Steps:
# 1 - Preprocess features
# 2 - Split data in train and test groups
#     taking into account the unbalanced data distribution
#     between bloom and normal categories in the outcome variable
# 3 - Define hyperparameters for model tuning.
# 4 - Train and tune a RF model with 5 k-fold CV repeated 100 times
#     exploring different mtry hyperparameters
#     using parallel computing.
# 5 - Train and tune RF models using 100 data splits and 
#     5 and 10 k-folds
#---------------------------------------------

############################################
# Load dependencies
############################################

library(tidyverse)
library(mikropml)
library(caret)

library(tictoc)
library(future)
library(doFuture)
library(future.apply)

set.seed(20210401)

load("processed_data/asvs_filtered.RData")

# we are working with proportion data
# remove other data to avoid mistakes
rm("asvs_co","asvs_co_wide")

############################################
# Preprocess features
############################################

# Preprocess features by removing near-zero variance 
# and highly correlated features

# randomize data row-wise (not really needed)
asvs_ra_wide <- asvs_ra_wide[sample(nrow(asvs_ra_wide)),]

asvs_prep <- asvs_ra_wide %>% 
  select(-sample, -gr_event) %>%
  mutate(st = as.factor(st),
         year = as.factor(year),
         month = as.factor(month),
         season = as.factor(season),
         depth_lev = as.factor(depth_lev)) %>%
  mikropml::preprocess_data(outcome_colname = "event",
                  method = c("center","scale"),
                  to_numeric = FALSE,
                  collapse_corr_feats = TRUE, 
                  group_neg_corr = TRUE,    
                  remove_var='nzv',          
                  prefilter_threshold = -1) # do not filter by prevalence

# check number of removed features
length(asvs_prep$removed_feats)
#[1] 184

# get preprocessed data 
asvs <- asvs_prep$dat_transformed

# Boxplot preprocessed data
pdf(file="results/boxplots/asvs_boxplot_preprocessed_mod1.pdf")
asvs %>%
  pivot_longer(-event, names_to = "feature", values_to = "value") %>%
  ggplot(aes(y=value, x=feature)) +
  geom_boxplot(outlier.size = 0.4) +
  theme(axis.text.x = element_blank())
dev.off()

# save list of removed features
write.csv2(asvs_prep$removed_feats, file = "results/asvs_features_removed_mod1.csv")

##########################################################
# Data split 
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data
train_part <- caret::createDataPartition(
  y=asvs_ra_wide$event,
  times = 100, 
  p = 0.8)

# Get the first split for the initial test
train_indices <- train_part[[1]]

# Check train data distribution among groups
lapply(train_part, function(x, y) table(y[x]), y = asvs_ra_wide$event)
# $Resample1
# 
# bloom normal 
# 25    108 


# Check test data distribution among groups
test_part <-list()
for(i in 1:length(train_part)){
  test_part[[i]] <-which(!(1:nrow(asvs_ra_wide) %in% train_part[[i]]))
}
lapply(test_part, function(x, y) table(y[x]), y = asvs_ra_wide$event)
# [[1]]
# 
# bloom normal 
# 6     27 


###################################################
# Model training - single split - inital test
##################################################

# check default hyperparameters
get_hyperparams_list(asvs,"rf")
# $mtry
# [1] 16 32 64


# Model train
# single data split with 5 kfold cross-validation X 100
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 10)

tic()
results_rf1 <- run_ml(asvs,
         method = "rf",
         outcome_colname = "event",
         training_frac = train_indices,
         kfold = 5,
         cv_times = 100,
         find_feature_importance = FALSE,
         seed = 20210401)
toc()

future::plan(future::sequential)

# Check mtry evaluation
performance <- get_hp_performance(results_rf1$trained_model)
plot_hp_performance(performance$dat, mtry, AUC)


# expand grid search for mtry
tune_grid <- list(mtry = c(8,16,32,48,64))

# train model again with expanded grid
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 10)

tic()
results_rf1 <- run_ml(asvs,
                      method = "rf",
                      outcome_colname = "event",
                      training_frac = train_indices,
                      kfold = 5,
                      cv_times = 100,
                      hyperparameters = tune_grid,
                      find_feature_importance = FALSE,
                      seed = 20210401)
toc()


# Model performance
performance <- get_hp_performance(results_rf1$trained_model)
plot_hp_performance(performance$dat, mtry, AUC)


calc_perf_metrics(results_rf1$test_data,
                  results_rf1$trained_model,
                  "event",
                  twoClassSummary,
                  class_probs = TRUE)

# ROC      Sens      Spec 
# 0.9135802 0.9629630 0.5000000

saveRDS(results_rf1, file = "results/models/asvs_rf1.RDS")


################################################
# Model training - 100 data splits
################################################

# 5 k-fold
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 100)

tic()
results_multi <- future.apply::future_lapply(train_part, function(train_ind){
  run_ml(asvs,
         method = "rf",
         outcome_colname = "event",
         training_frac = train_ind,
         kfold = 5,
         cv_times = 10,
         hyperparameters = tune_grid,
         find_feature_importance = FALSE,
         seed = 20210401)
  }, future.seed = TRUE)
toc()

saveRDS(results_multi, file = "results/models/asvs_rf1_multi_k5.RDS")

# 10 k-fold
tic()
results_multi_2 <- future.apply::future_lapply(train_part, function(train_ind){
  run_ml(asvs,
         method = "rf",
         outcome_colname = "event",
         training_frac = train_ind,
         kfold = 10,
         cv_times = 10,
         hyperparameters = tune_grid,
         find_feature_importance = FALSE,
         seed = 20210401)
  }, future.seed = TRUE)
toc()

saveRDS(results_multi_2, file = "results/models/asvs_rf1_multi_k10.RDS")

# compare performances
perf_dfk5 <- future.apply::future_lapply(results_multi,
  function(result) {
    result[["performance"]] %>%
      select(cv_metric_AUC, AUC,F1, method)
  },
  future.seed = TRUE
) %>%
  dplyr::bind_rows()

perf_dfk10 <- future.apply::future_lapply(results_multi_2,
  function(result) {
    result[["performance"]] %>%
      select(cv_metric_AUC, AUC,F1, method)
  },
  future.seed = TRUE
) %>%
  dplyr::bind_rows()



future::plan(sequential)
 


