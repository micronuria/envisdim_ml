############################################
# Load dependencies
############################################

library(tidyverse)
library(mikropml)
library(caret)
# 
# library(tictoc)
# library(future)
# library(doFuture)
# library(future.apply)
# 
# set.seed(20210401)

# models with several data splits between train and test

multi <- readRDS("results/models/clusters_rf1_multi_k10.RDS")

# model <- cl_rf1_k10[[5]]$trained_model
# test_data <- cl_rf1_k10[[5]]$test_data
# confusionMatrix(predict(model, test_data),as.factor(test_data$event), positive = "bloom")
# #####
# 
# 
# models <- lapply(cl_rf1_k10, function(x) x$trained_model)
# hp_metrics <- combine_hp_performance(models)
# plot_hp_performance(hp_metrics$dat, mtry, AUC)



lapply(multi, function(x) { x$trained_model$results }) %>% 
  dplyr::bind_rows() %>%
  select(mtry, Accuracy, AUC, F1, Kappa, Precision, Recall, Sensitivity, Specificity) %>%
  pivot_longer(-mtry, values_to = "value", names_to = "index") %>%
   ggplot(aes(x=mtry, y=value, group = mtry)) +
   ggtitle("cl_ra_multi_k10") +
     geom_boxplot() +
  facet_wrap(~index)
