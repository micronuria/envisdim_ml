###########################
#
# MODEL A TRAINING
# clusters + proportion
#
###########################

library(tidyverse)
library(mikropml)
#library(purr)
library(future)
library(furrr)

library(caret)

library(tictoc)
set.seed(20210401)

load("processed_data/clusters_preprocess.RData")


# Data split
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data

train_part <- caret::createDataPartition(
  y=composite_wide$gr_event,
  times = 1, 
  p = 0.8)

# Check train data distribution among groups
lapply(train_part, function(x, y) table(y[x]), y = composite_wide$gr_event)

# Check test data distribution among groups
test_part <-list()
for(i in 1:length(train_part)){
  test_part[[i]] <-which(!(1:nrow(composite_wide) %in% train_part[[i]]))
}
lapply(test_part, function(x, y) table(y[x]), y = composite_wide$gr_event)

# split transfomed data
clusters <- clusters_prep$dat_transformed
train_data <- clusters[train_part[[1]],]
test_data <- clusters[-train_part[[1]],]

###################################################
# for k-fold splits, use the standard random split 
# within the event variable (normal-bloom) 
# and ignore the campaigs - number of samples is limited.

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
prueba <- run_ml(train_data,
       method = "rf",
       outcome_colname = "event",
       kfold=5,
       cv_times = 100,
       seed = 1)
toc()




tic()
iterative_results <- furrr::future_map(c(1), get_results)
toc()

future::plan(sequential)



