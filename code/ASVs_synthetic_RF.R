
# Author: Nuria Fernández-Gonzalez

###########################################
#
# Synthetic ASVs with SMOTE
#
###########################################

#---------------------------------------------
# Description. This script trains RF models
# using synthetic data generation

# Steps:
# 1 - Preprocess ASVs data
# 2 - Define divisions of data in train and test groups
#     taking into account the unbalanced data distribution
#     between bloom and normal categories in the outcome variable
# 3 - Define hyperparameters for model tuning.
# 4 - Define train details: tune a RF model with 5 k-fold CV repeated 10 times
#     and generate synthetic data with SMOTE
# 5 - Train and tune RF models using 100 data splits 
#     using parallel computing
#---------------------------------------------

############################################
# Load dependencies
############################################

library(tidyverse)
library(mikropml)
library(caret)
library(themis)
library(tictoc)
library(future)
library(doFuture)
library(future.apply)


set.seed(20210401)

load("processed_data/asvs_filtered.RData")

# we are working with proportion data
rm("asvs_co","asvs_co_wide")

############################################
# Preprocess ASVs
############################################

# Preprocess ASVs by removing near-zero variance 
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
                  prefilter_threshold = -1) 

# get preprocessed data 
asvs <- asvs_prep$dat_transformed
asvs$event <- factor(asvs$event, levels = c("bloom", "normal"))

# 100 data splits
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data
train_list <- caret::createDataPartition(
  y=asvs$event,
  times = 100,
  p = 0.80)


# Define hyperparameters grid and Train details
##########################################################

# hyperparameters grid search
tune_grid <- expand.grid(mtry = c(8,16,32,48,64))

# Define 5 k-fold RCV 10 times
# generate synthetic data with SMOTE algorithm
train_control <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 10,
  sampling = "smote",
  classProbs = TRUE,
  verboseIter = FALSE)

# Train models for 100 data splits
############################################

# start results object
results <- list()

# paralellize computing
doFuture::registerDoFuture()
future::plan(future::multicore, workers = 10)

tic()
for (i in seq_len(length(train_list))) {
  # generate train and test groups
  train_indices <- train_list[[i]]
  train_group <- asvs[train_indices,]
  test_indices <- !(1:nrow(asvs) %in% train_list[[i]])
  test_group <- asvs[test_indices,]
  
  # train model 
  modelo <- train(event ~ .,
                     data = train_group,
                     method = "rf",
                     trControl = train_control,
                     tuneGrid = tune_grid,
                     metric = "ROC",
                     trace = FALSE)
  
  # Make predicitions with test data
  pred_class <- predict(modelo_rf_1, newdata = test_group)
  
  # Calculate model performance
  performance <- confusionMatrix(data=pred_class, test_group$event,
                            positive = "bloom", mode = "everything")
  
  # join results
  results[[i]] <- list(
    trained_model = modelo,
    test_data = test_group,
    performance = performance
    )
}
toc()

future::plan(future::sequential)

# Save results
saveRDS(results, "results/models/test_syntetic_results.RDS")


