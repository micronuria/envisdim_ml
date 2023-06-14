

# Author: Nuria Fernández-Gonzalez

###########################################
#
#  RF model with syntetic ASVs
#
###########################################

#---------------------------------------------
# Description. This script trains a RF model
# using 
# - ASVs and some environmental variables

#---------------------------------------------


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
asvs$event <- factor(asvs$event, levels = c("bloom", "normal"))

# Data split
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data
train_list <- caret::createDataPartition(
  y=asvs$event,
  times = 1, 
  p = 0.80)

# split transfomed data
train_indices <- train_list[[1]]
train_group <- asvs[train_indices,]
test_indices <- !(1:nrow(asvs) %in% train_list[[1]])
test_group <- asvs[test_indices,]


# Check data distribution among event clases
table(train_group$event)
# bloom normal 
# 25    108 

table(test_group$event)
# bloom normal 
# 6     27 

# Create synthetic data
train.themis <- themis::smote(train_group, var = "event", over_ratio = 3)
table(train.themis$event)
# bloom normal 
#   324    324 

# Model train
# hyperparameters grid search
tune_grid <- expand.grid(mtry = c(8,16,32,48,64))

# train control
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 80,
  verboseIter = FALSE)

# Model train
# single data split with 10 k-fold cross-validation X 10
doFuture::registerDoFuture()
future::plan(future::multicore, workers = 80)

tic()
modelo_rf_1 <- train(event ~ .,
                     data = train.themis,
                     method = "rf",
                     trControl = train_control,
                     tuneGrid = tune_grid,
                     metric = "Kappa",
                     trace = FALSE)
                   
toc()

future::plan(future::sequential)

saveRDS(modelo_rf_1, file = "results/models/asvs_rf1_smote.RDS")


# Check model performance
# Class predictions
pred_class <- predict(modelo_rf_1, newdata = test_group)
perf_rf1 <- confusionMatrix(data=pred_class, test_group$event,
                            positive = "bloom", mode = "everything")
perf_rf1

