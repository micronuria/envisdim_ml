
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# ASVs RF MODEL 2
#
###########################################

#---------------------------------------------
# Description. This script trains a RF model
# using 
# - ASVs and some environmental variables.
# - ASVs NOT filtered by abundance.
# - Normalization method: proportions.
# - Not filtered by prevalence.
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

load(file="processed_data/asvs_no_filtered.RData")


# randomize data row-wise (not really needed)
asvs_ra_wide <- asvs_ra_wide[sample(nrow(asvs_ra_wide)),]

asvs_prep <- asvs_ra_wide %>% 
  select(-sample, -gr_event) %>%
  mikropml::preprocess_data(outcome_colname = "event",
                  method = c("center","scale"),
                  collapse_corr_feats = TRUE, 
                  group_neg_corr = FALSE,    
                  remove_var='nzv',          
                  prefilter_threshold = -1) # do not filter by prevalence

length(asvs_prep$removed_feats)
#[1] 3515

# get preprocessed data 
asvs <- asvs_prep$dat_transformed

# Boxplot preprocessed data
pdf(file="results/boxplots/asvs_boxplot_preprocessed_2.pdf")
asvs %>%
  pivot_longer(-event, names_to = "feature", values_to = "datos") %>%
  ggplot(aes(y=datos, x=feature)) +
  geom_boxplot(outlier.size = 0.4) +
  theme(axis.text.x = element_blank())
dev.off()

##########################################################
# Data split and hyperparameters definition
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data
train_part <- caret::createDataPartition(
  y=asvs_ra_wide$event,
  times = 100, 
  p = 0.8)

# split transfomed data
train_indices <- train_part[[1]]

# define hyperparameters to explore
get_hyperparams_list(asvs,"rf")
# $mtry
# [1] 32 64 128

tune_grid <- list(mtry = c(16,32,48,64,128))

###################################################
# Model training - single split
##################################################
# Model train
# single data split with 5 cross-validation X 100

doFuture::registerDoFuture()
future::plan(future::multicore, workers = 100)

tic()
results_rf2 <- run_ml(asvs,
         method = "rf",
         outcome_colname = "event",
         training_frac = train_indices,
         kfold = 5,
         cv_times = 100,
         hyperparameters = tune_grid,
         find_feature_importance = FALSE,
         seed = 20210401)
toc()

saveRDS(results_rf2, file = "results/models/asvs_rf2_k5.RDS")

################################################
# Model training - 100 data splits
################################################

# 5 k-fold
tic()
results_multi_k5 <- future.apply::future_lapply(train_part, function(train_ind){
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

saveRDS(results_multi_k5, file = "results/models/asvs_rf2_multi_k5.RDS")

# 10 k-fold x10 
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 100)

tic()
results_multi_k10 <- future.apply::future_lapply(train_part, function(train_ind){
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

saveRDS(results_multi_k10, file = "results/models/asvs_rf2_multi_k10.RDS")


save(list=ls(), file="processed_data/asvs_model_2.RData")
 