

############################################
# Load dependencies
############################################

library(tidyverse)
library(caret)
library(mikropml)
#library(ROSE)
library(themis)

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
  mutate(event = as.factor(event)) %>%
  #rename_with(., ~ gsub("-", "_", .x, fixed = TRUE)) %>%
  mikropml::preprocess_data(outcome_colname = "event",
                            method = c("center","scale"),
                            collapse_corr_feats = TRUE, 
                            group_neg_corr = FALSE,    
                            remove_var='nzv',          
                            prefilter_threshold = -1) # do not filter by prevalence

length(clusters_prep$removed_feats)
#[1] 394

# get preprocessed data and randomize (not really needed, caret randomizes the split)
clusters <- clusters_prep$dat_transformed

##########################################################
# Data splits - synthetic data
##########################################################

# check class balance
table(clusters$event)
# bloom normal 
# 31    135

# Create a split of data ensuring samples from 
# both events and campaigns are present in train and test data
train_list <- caret::createDataPartition(
  y=clusters_ra_wide$event,
  times = 100, 
  p = 0.8)

# split transfomed data
train_indices <- train_list[[1]]
train_group <- clusters[train_indices,]
test_indices <- !(1:nrow(clusters) %in% train_list[[1]])
test_group <- as.data.frame(clusters[test_indices,])

# Check data distribution among event clases
table(train_group$event)
# bloom normal 
# 25    108 

table(test_group$event)
# bloom normal 
# 6     27 

# Create synthetic data
train.themis <- themis::smote(train_group, var = "event", over_ratio = 1)
table(train.themis$event)

###################################################
# Model training - single split
##################################################

# define hyperparameters to explore
# sqrt(p) ~ 18  c(3,6,9,18,24,36)
tune_grid <- expand.grid(mtry = c(3,6,9,18,24,36))

# train control
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  verboseIter = FALSE)

# Model train
# single data split with 10 k-fold cross-validation X 10
doFuture::registerDoFuture()
future::plan(future::multisession, workers = 10)

tic()
modelo_rf_1 <- train(event ~ .,
                     data = train.themis,
                     method = "rf",
                     trControl = train_control,
                     tuneGrid = tune_grid,
                     metric = "Kappa",
                     trace = FALSE)
                   

toc()

saveRDS(results, file = "results/models/clusters_rf1_smote.RDS")

future::plan(future::sequential)


# Check model performance

# Class predictions
pred_class <- predict(modelo_rf_1, newdata = test_group)
perf_rf1 <- confusionMatrix(data=pred_class, test_group$event,
                            positive = "bloom", mode = "everything")
perf_rf1
