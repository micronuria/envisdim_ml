
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

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
# and highly correlated features, but not negatively correlated features
# for interpretability

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
#[1] 184

# get preprocessed data 
asvs <- asvs_prep$dat_transformed

# Boxplot preprocessed data
pdf(file="results/boxplots/asvs_boxplot_preprocessed_1.pdf")
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

# define hyperparameters to explore
get_hyperparams_list(asvs,"rf")
# $mtry
# [1] 16 32 64

tune_grid <- list(mtry = c(8,16,32,48,64))

###################################################
# Model training - single split
##################################################


# Model train
# single data split with 5 cross-validation X 100

doFuture::registerDoFuture()
future::plan(future::multicore, workers = 100)

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

saveRDS(results_rf1, file = "results/models/asvs_rf1.RDS")

future::plan(future::sequential)

# Model performance
performance <- get_hp_performance(results_rf1$trained_model)
plot_hp_performance(performance$dat, mtry, AUC)

results_rf1$trained_model
#Random Forest 
#
# 133 samples
#1031 predictors
#   2 classes: 'bloom', 'normal' 
#
#No pre-processing
#Resampling: Cross-Validated (5 fold, repeated 100 times) 
#Summary of sample sizes: 106, 107, 106, 107, 106, 106, ... 
#Resampling results across tuning parameters:
#
#  mtry  logLoss    AUC        prAUC      Accuracy   Kappa      F1       
#   8    0.3446908  0.8742970  0.7109381  0.8514017  0.3668316  0.4826984
#  16    0.3425114  0.8743238  0.7107516  0.8489088  0.3634849  0.4816698
#  32    0.3414061  0.8737364  0.7082347  0.8455214  0.3621877  0.4797427
#  48    0.3411306  0.8729970  0.7076600  0.8439972  0.3613582  0.4853958
#  64    0.3412128  0.8722926  0.7056209  0.8440085  0.3645412  0.4884509
#  Sensitivity  Specificity  Pos_Pred_Value  Neg_Pred_Value  Precision  Recall
#  0.3300       0.9721472    0.7726190       0.8638021       0.7726190  0.3300
#  0.3340       0.9681515    0.7550529       0.8640720       0.7550529  0.3340
#  0.3432       0.9618442    0.7260923       0.8649617       0.7260923  0.3432
#  0.3496       0.9584805    0.7059834       0.8658398       0.7059834  0.3496
#  0.3548       0.9572944    0.7038733       0.8666420       0.7038733  0.3548
#  Detection_Rate  Balanced_Accuracy
#  0.06205413      0.6510736        
#  0.06280342      0.6510758        
#  0.06453561      0.6525221        
#  0.06573789      0.6540403        
#  0.06671510      0.6560472        
#
#
#AUC was used to select the optimal model using the largest value.
#The final value used for the model was mtry = 16.

calc_perf_metrics(results_rf1$test_data,
                  results_rf1$trained_model,
                  "event",
                  twoClassSummary,
                  class_probs = TRUE)

#      ROC      Sens      Spec 
#0.9043210 0.9629630 0.3333333 


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

mean(perf_dfk5$F1)
#[1] 0.883613

mean(perf_dfk10$F1)
#[1] 0.885509
mean(perf_dfk5$AUC)
#[1] 0.8695988
mean(perf_dfk10$AUC)
#[1] 0.8714506

future::plan(sequential)
 

save(list=ls(), file="processed_data/asvs_model_1.RData")
