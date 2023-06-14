# Author: Nuria Fernández-Gonzalez

###########################################
#
# ASVs- SVM classification
#
###########################################

#---------------------------------------------
# Description. This script trains a SVM model
# using
# - ASVs and some environmental variables
# - Normalization method: proportions
# - Not filtered by prevalence
# Steps:
# 1 - Preprocess asvs data
# 2 - Split data in train and test groups
#     taking into account the unbalanced data distribution
#     between bloom and normal categories in the outcome variable
# 3 - Train and tune a SVMRadial model with 5 k-fold CV repeated 10 times
#     exploring different hyperparameters
#     using parallel computing. 100 data splits
# 4 - Obtain cross-validation results
#---------------------------------------------

##########################################
# Load dependencies
############################################

library(tidyverse)
library(caret)

library(mikropml)

library(tictoc)
library(future)
library(doFuture)
library(future.apply)

library(ggpubr)

set.seed(20210401)

load("processed_data/asvs_filtered.RData")

# we are working with proportion data
rm("asvs_co","asvs_co_wide")

############################################
# Preprocess ASVs
############################################

# Preprocess ASVs by removing near-zero variance 
# and highly correlated features, but not negatively correlated features
# for interpretability

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

length(asvs_prep$removed_feats)
#[1] 184

# get preprocessed data 
asvs <- asvs_prep$dat_transformed

##########################################################
# Data split and define hyperparameters grid
#########################################################
# Create balanced split of data - ensure samples from
# both events and campaigns are present in train and test data
train_part <- caret::createDataPartition(
  y=asvs_ra_wide$event,
  times = 100, 
  p = 0.8)

# define hyperparameters to explore
# Using defaults:
get_hyperparams_list(asvs, "svmRadial")
#$C
#[1] 1e-03 1e-02 1e-01 1e+00 1e+01 1e+02

#$sigma
#[1] 1e-06 1e-05 1e-04 1e-03 1e-02 1e-01


################################################
# Model training - 100 data splits
################################################

doFuture::registerDoFuture()
future::plan(future::multicore, workers = 100)

tic()
asvs_svm1_k5 <- future.apply::future_lapply(train_part, function(train_ind){
  run_ml(asvs,
         method = "svmRadial",
         outcome_colname = "event",
         training_frac = train_ind,
         kfold = 5,
         cv_times = 10,
         find_feature_importance = FALSE,
         seed = 20210401)
  }, future.seed = TRUE)
toc()

saveRDS(asvs_svm1_k5, file = "results/models/asvs_svm1_multi100_k5_10.RDS")


#########################################################
# 
asvs_svm1_k5 <- readRDS("results/models/asvs_svm1_multi100_k5_10.RDS")


# Validation results

models <- lapply(asvs_svm1_k5, function(x) x$trained_model)
hp_metrics <- combine_hp_performance(models)
pc <- plot_hp_performance(hp_metrics$dat, C, AUC)
pc <- pc + scale_x_continuous(trans='log10') 
ps <- plot_hp_performance(hp_metrics$dat, sigma, AUC)
ps <- ps + scale_x_continuous(trans='log10') 

ggarrange(pc, ps,
          labels = c("a)", "b)"),
          ncol = 1, nrow = 2)


# Performance results
perf_df <- lapply(asvs_svm1_k5, function(result) {
                  result[["performance"]] %>%
    select(cv_metric_AUC, AUC, Accuracy, Kappa, Sensitivity, Specificity)}) %>%
  dplyr::bind_rows() 

summary(perf_df)
# 
# cv_metric_AUC         AUC            Accuracy          Kappa          Sensitivity    
# Min.   :0.8077   Min.   :0.5123   Min.   :0.7273   Min.   :-0.1000   Min.   :0.0000  
# 1st Qu.:0.8502   1st Qu.:0.8194   1st Qu.:0.8182   1st Qu.: 0.0000   1st Qu.:0.8889  
# Median :0.8642   Median :0.8765   Median :0.8182   Median : 0.2466   Median :0.9630  
# Mean   :0.8649   Mean   :0.8597   Mean   :0.8339   Mean   : 0.2428   Mean   :0.9026  
# 3rd Qu.:0.8810   3rd Qu.:0.9275   3rd Qu.:0.8485   3rd Qu.: 0.4563   3rd Qu.:1.0000  
# Max.   :0.9178   Max.   :0.9877   Max.   :0.9091   Max.   : 0.7130   Max.   :1.0000  
# Specificity    
# Min.   :0.0000  
# 1st Qu.:0.0000  
# Median :0.2500  
# Mean   :0.3228  
# 3rd Qu.:0.5417  
# Max.   :1.0000  


# boxplot
perf_df %>% 
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) 




