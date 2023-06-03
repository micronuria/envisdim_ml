
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# RF using clusters with rclr
# Model 2
#
###########################################

#---------------------------------------------
# Description. This script trains a RF model
# using 
# - clusters and some environmental variables
# - Normalization method: rclr
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
rm("clusters_ra", "clusters_ra_wide")

############################################
# Preprocess clusters
############################################

# Preprocess genus clusters by removing near-zero variance 
# and highly correlated features, but not negatively correlated features
# for interpretability
# Also removing features not present in at least two samples

library(vegan)

# prepare data
# apply robust center log ratio transformation to count data
clusters_rclr <- decostand(clusters_co_wide[,11:ncol(clusters_co_wide)],
                         "rclr", MARGIN = 1)

clusters_co_wide[,11:ncol(clusters_co_wide)] <- clusters_rclr

# randomize data row-wise (not really needed)
clusters_co_wide <- clusters_co_wide[sample(nrow(clusters_co_wide)),]

clusters_prep <- clusters_co_wide %>% 
  select(-sample, -gr_event) %>%
  mikropml::preprocess_data(outcome_colname = "event",
                            method = c("center","scale"),
                            collapse_corr_feats = TRUE, 
                            group_neg_corr = FALSE,    
                            remove_var='nzv',          
                            prefilter_threshold = -1) 

length(clusters_prep$removed_feats)
#[1] 392

# get preprocessed data and randomize (not really needed, caret randomizes the split)
clusters <- clusters_prep$dat_transformed

# Boxplot preprocessed data
pdf(file="results/boxplots/clusters_boxplot_preprocessed_2.pdf")
clusters %>%
  pivot_longer(-event, names_to = "feature", values_to = "datos") %>%
  ggplot(aes(y=datos, x=feature)) +
  geom_boxplot(outlier.size = 0.4) +
  theme(axis.text.x = element_blank())
dev.off()

##########################################################
# Data split - hyperparameters setting
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data
train_part <- caret::createDataPartition(
  y=clusters_co_wide$event,
  times = 100, 
  p = 0.8)

# split transfomed data
train_indices <- train_part[[1]]

# Check train data distribution among groups
lapply(train_part, function(x, y) table(y[x]), y = clusters_co_wide$event)
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


# define hyperparameters to explore
# sqrt(p) ~ 18
tune_grid <- list(mtry = c(3,6,9,18,24,36))

###################################################
# Model training - single split
##################################################
# Model train
# single data split with 5 cross-validation X 100

doFuture::registerDoFuture()
future::plan(future::multisession, workers = 100)

tic()
results_rf2 <- run_ml(clusters,
                      method = "rf",
                      outcome_colname = "event",
                      training_frac = train_indices,
                      kfold = 5,
                      cv_times = 100,
                      hyperparameters = tune_grid,
                      find_feature_importance = FALSE,
                      seed = 20210401)
toc()

saveRDS(results_rf2, file = "results/models/clusters_rf2.RDS")

future::plan(future::sequential)


# Model performance
performance <- get_hp_performance(results_rf2$trained_model)
plot_hp_performance(performance$dat, mtry, AUC)

calc_perf_metrics(results_rf2$test_data,
                  results_rf2$trained_model,
                  "event",
                  twoClassSummary,
                  class_probs = TRUE)


################################################
# Model training - 100 data splits
################################################

# 5-kfold
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 100)

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

saveRDS(results_multi2, file = "results/models/clusters_rf2_multi_k5.RDS")

# 10 k-folds
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 100)

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

saveRDS(results_multi_k10, file = "results/models/clusters_rf2_multi_k10.RDS")

future::plan(sequential)

###############################
# save all data
save(list=ls(), file="processed_data/clusters_model_2.RData")
