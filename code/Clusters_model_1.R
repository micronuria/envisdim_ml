
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# CLUSTERS PREPROCESS
# for model 1
#
###########################################

#---------------------------------------------
# Description. This script:
# 1 - Preprocess clusters data
#
# 2 - 
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


############################################
# Preprocess clusters
############################################

# Preprocess genus clusters by removing near-zero variance 
# and highly correlated features, but not negatively correlated features
# for interpretability

clusters_prep <- clusters_ra_wide %>% 
  select(-sample, -gr_event) %>%
  preprocess_data(outcome_colname = "event",
                  method = c("center","scale"),
                  collapse_corr_feats = TRUE, 
                  group_neg_corr = FALSE,    
                  remove_var='nzv',          
                  prefilter_threshold = -1) # do not filter by prevalence




# Data split
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data

train_part <- caret::createDataPartition(
  y=clusters_ra_wide$gr_event,
  times = 1, 
  p = 0.8)

# Check train data distribution among groups
lapply(train_part, function(x, y) table(y[x]), y = clusters_ra_wide$gr_event)
# $Resample1
# 
# bloom_dim  bloom_env normal_dim normal_env 
# 8         17         21         88 



# Check test data distribution among groups
test_part <-list()
for(i in 1:length(train_part)){
  test_part[[i]] <-which(!(1:nrow(clusters_ra_wide) %in% train_part[[i]]))
}
lapply(test_part, function(x, y) table(y[x]), y = clusters_ra_wide$gr_event)
# [[1]]
# 
# bloom_dim  bloom_env normal_dim normal_env 
# 2          4          5         21 

# split transfomed data
clusters <- clusters_prep$dat_transformed
train_data <- clusters[train_part[[1]],]
test_data <- clusters[-train_part[[1]],]

###################################################
# for k-fold splits, use the standard random split 
# within the event variable (normal-bloom) 
# and ignore the campaigs - number of samples is limited.

### USAR CARET



# Model train


# Define function with model
get_results <- function(seed) {
  run_ml(train_data,
         method = "rf",
         outcome_colname = "event",
         kfold=5,
         cv_times = 10,
         seed = seed,
  )
}

# check time for a single run with 5 cross-validation X 100

doFuture::registerDoFuture()
future::plan(future::multisession, workers = 10)

tic()
results_rf1 <- run_ml(train_data,
                 method = "rf",
                 outcome_colname = "event",
                 kfold=5,
                 cv_times = 100,
                 seed = 1)
toc()

future::plan(future::sequential)


tic()
iterative_results <- furrr::future_map(c(1), get_results)
toc()

future::plan(sequential)




