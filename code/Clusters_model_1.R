
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# RF using clusters with proportions
# Model 1
#
###########################################

#---------------------------------------------
# Description. This script trains RF models
# using 
# - clusters and some environmental variables
# - Normalization method: proportions
# - Not filtered by prevalence
# Steps:
# 1 - Preprocess clusters data
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

load("processed_data/clusters_grouping.RData")

# we are working with proportion data
# remove other data to avoid mistakes
rm("clusters_count", "clusters_co_wide")

############################################
# Preprocess clusters
############################################

# Preprocess genus clusters by removing near-zero variance 
# and highly correlated features

# randomize data row-wise (not really needed)
clusters_ra_wide <- clusters_ra_wide[sample(nrow(clusters_ra_wide)),]

clusters_prep <- clusters_ra_wide %>% 
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
                  remove_var= 'nzv',          
                  prefilter_threshold = -1) # do not filter by prevalence

length(clusters_prep$removed_feats)
#[1] 394

# get preprocessed data and randomize (not really needed, caret randomizes the split)
clusters <- clusters_prep$dat_transformed

# Boxplot preprocessed data
pdf(file="results/boxplots/clusters_boxplot_preprocessed_mod1.pdf")
clusters %>%
  pivot_longer(-event, names_to = "feature", values_to = "value") %>%
  ggplot(aes(y=value, x=feature)) +
  geom_boxplot(outlier.size = 0.4) +
  theme(axis.text.x = element_blank())
dev.off()

# save list of removed features
write.csv2(clusters_prep$removed_feats, file = "results/clusters_features_removed_mod1.csv")

##########################################################
# Data splits - define hyperparameters
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data
train_part <- caret::createDataPartition(
  y=clusters_ra_wide$event,
  times = 100, 
  p = 0.8)

# split transfomed data
train_indices <- train_part[[1]]

# Check train data distribution among groups
lapply(train_part, function(x, y) table(y[x]), y = clusters_ra_wide$event)
# $Resample1
# 
# bloom normal 
# 25    108 


# Check test data distribution among groups
test_part <-list()
for(i in 1:length(train_part)){
  test_part[[i]] <-which(!(1:nrow(clusters_ra_wide) %in% train_part[[i]]))
}
lapply(test_part, function(x, y) table(y[x]), y = clusters_ra_wide$event)
# [[1]]
# 
# bloom normal 
# 6     27 


###################################################
# Model training - single split
##################################################

# define hyperparameters to explore
get_hyperparams_list(clusters, "rf")
#$mtry
# [1]  9 18 36

# Model train
# single data split with 5 cross-validation X 100
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 10)

tic()
results_rf1 <- run_ml(clusters,
         method = "rf",
         outcome_colname = "event",
         training_frac = train_indices,
         kfold = 5,
         cv_times = 100,
         find_feature_importance = FALSE,
         seed = 20210401)
toc()

future::plan(future::sequential)

# Check model performance to adjust mtry grid search
performance <- get_hp_performance(results_rf1$trained_model)
plot_hp_performance(performance$dat, mtry, AUC)

# Adjust mtry searh and repeat training with a single split
tune_grid <- list(mtry = c(3,6,9,18,24,36))

doFuture::registerDoFuture()
future::plan(future::multisession, workers = 10)

tic()
results_rf1 <- run_ml(clusters,
                      method = "rf",
                      outcome_colname = "event",
                      training_frac = train_indices,
                      kfold = 5,
                      cv_times = 100,
                      hyperparameters = tune_grid,
                      find_feature_importance = FALSE,
                      seed = 20210401)
toc()

future::plan(future::sequential)

# Check model performance again
performance <- get_hp_performance(results_rf1$trained_model)
plot_hp_performance(performance$dat, mtry, AUC)


saveRDS(results_rf1, file = "results/models/clusters_single_rf1_k5.RDS")


################################################
# Model training - 100 data splits
################################################

doFuture::registerDoFuture()
future::plan(future::multisession, workers = 100)
# 5 k-folds
tic()
results_multi_k5 <- future.apply::future_lapply(train_part, function(train_ind){
  run_ml(clusters,
         method = "rf",
         outcome_colname = "event",
         training_frac = train_ind,
         kfold = 5,
         cv_times = 100,
         hyperparameters = tune_grid,
         find_feature_importance = FALSE,
         seed = 20210401)
  }, future.seed = TRUE)
toc()

saveRDS(results_multi_k5, file = "results/models/clusters_rf1_multi100_k5_100.RDS")

# 10 k-folds
tic()
results_multi_k10 <- future.apply::future_lapply(train_part, function(train_ind){
  run_ml(clusters,
         method = "rf",
         outcome_colname = "event",
         training_frac = train_ind,
         kfold = 10,
         cv_times = 100,
         hyperparameters = tune_grid,
         find_feature_importance = FALSE,
         seed = 20210401)
  }, future.seed = TRUE)
toc()

future::plan(sequential)

saveRDS(results_multi_k10, file = "results/models/clusters_rf1_multi100_k10_100.RDS")
