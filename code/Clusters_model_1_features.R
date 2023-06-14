
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# CLUSTERS model 1
# features
###########################################

#---------------------------------------------
# Description. This script trains a RF model
# using 

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
# and highly correlated features, but not negatively correlated features
# for interpretability

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

# Data split
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data
train_part <- caret::createDataPartition(
  y=clusters_ra_wide$event,
  times = 100, 
  p = 0.8)

# define hyperparameters to explore
tune_grid <- list(mtry = c(6))

################################################
# Model training - 100 data splits
################################################

doFuture::registerDoFuture()
future::plan(future::multisession, workers = 100)

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
         find_feature_importance = TRUE,
         seed = 20210401)
  }, future.seed = TRUE)
toc()

saveRDS(results_multi_k10, file = "results/models/clusters_rf1_multi100_k10_100_features.RDS")

future::plan(sequential)

