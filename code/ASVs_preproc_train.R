

library(tidyverse)
library(mikropml)
#library(purr)
library(future)
library(furrr)

library(caret)

library(tictoc)

set.seed(20210401)

load("processed_data/initial_preprocess.RData")

# join count, taxonomy, metadata and ...
comp_asvs <- group_by(sample) %>%                                  # calculate relative abundances
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         # get rid of grouping structure on data
  select(-count) %>%                                    # get rid of count column
  inner_join(., chlclass, by="sample") %>%              # add data about chlorophyll groups
  inner_join(., met_raw)

# to wide format
comp_asvs_wide <- composite %>%
  pivot_wider(names_from = asv, values_from = rel_abun)


# Data split
##########################################################
# Create balanced split of data - ensure samples from 
# both events and campaigns are present in train and test data

train_part <- caret::createDataPartition(
  y=comp_asvs_wide$gr_event,
  times = 100, 
  p = 0.8)


# train
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



