
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# Initial model comparisons
#
###########################################

#---------------------------------------------
# Description: this script gathers together
# the performace of the different models
# developed using 100 data splits
# and performs an initial comparison
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

# models with several data splits between train and test
cl_rf1_k5 <- readRDS("results/models/clusters_rf1_multi_k5.RDS")
cl_rf1_k10 <- readRDS("results/models/clusters_rf1_multi_k10.RDS")

cl_rf2_k5 <- readRDS("results/models/clusters_rf2_multi_k5.RDS")
cl_rf2_k10 <- readRDS("results/models/clusters_rf2_multi_k10.RDS")

asvs_rf1_k5 <- readRDS("results/models/asvs_rf1_multi_k5.RDS")
asvs_rf1_k10 <- readRDS("results/models/asvs_rf1_multi_k10.RDS")

asvs_rf2_k5 <- readRDS("results/models/asvs_rf2_multi_k5.RDS")
asvs_rf2_k10 <- readRDS("results/models/asvs_rf2_multi_k10.RDS")


# get performance metrics for each set of models
get_performance <-function(modelo, x) {
  lapply(x, function(result) {result[["performance"]] %>% 
  select(cv_metric_AUC, AUC, Accuracy, Sensitivity, Specificity)}) %>%
  dplyr::bind_rows() %>%
  mutate(model = modelo)
}



perf_multi <- list()

perf_multi[[1]] <- get_performance("cl1-k5",cl_rf1_k5)
perf_multi[[2]] <- get_performance("cl1-k10",cl_rf1_k10) 
perf_multi[[3]] <- get_performance("cl2-k5",cl_rf2_k5)
perf_multi[[4]] <- get_performance("cl2-k10",cl_rf2_k10)
perf_multi[[5]] <- get_performance("asv1-k5",asvs_rf1_k5)
perf_multi[[6]] <- get_performance("asv1-k10",asvs_rf1_k10)
perf_multi[[7]] <- get_performance("asv2-k5",asvs_rf1_k5)
perf_multi[[8]] <- get_performance("asv2-k10",asvs_rf1_k10)

# make a boxplot to compare
Reduce(full_join,perf_multi) %>%
  pivot_longer(!model, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=model)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) #+
  #facet_wrap(~model)


###################
model <- cl_rf1_k10$Resample006

model$trained_model$results


###################
models <- lapply(cl_rf1_k10, function(x) x$trained_model)
hp_metrics <- combine_hp_performance(models)
plot_hp_performance(hp_metrics$dat, mtry, AUC)
##############################


# get performance metrics for each set of models
get_perf_mtry <-function(modelo, x) {
  lapply(x, function(result) {result$trained_model$results %>% 
      select(mtry, AUC, Accuracy, Kappa, Sensitivity, Specificity)}) %>%
    dplyr::bind_rows() %>%
    mutate(model = modelo)
}

df_perf <- get_perf_mtry("cl1-k10",asvs_rf1_k10)

plot_hp_performance(df_perf, mtry, Kappa)
plot_hp_performance(df_perf, mtry, Sensitivity)
