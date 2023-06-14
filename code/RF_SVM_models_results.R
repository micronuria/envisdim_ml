
# Author: Nuria Fernández-Gonzalez

###########################################
#
# RF and SVM model comparisons
#
###########################################

#---------------------------------------------
# Description: this script gathers together
# the validation and performace metrics
# of the RF and SVM models
# developed using 100 data splits
#---------------------------------------------


############################################
# Load dependencies
############################################

library(tidyverse)
library(mikropml)
library(ggpubr)
set.seed(20210401)

# RF models with several data splits between train and test
# clusters
cl_rf1_k5 <- readRDS("results/models/clusters_rf1_multi100_k5_100.RDS")
cl_rf1_k10 <- readRDS("results/models/clusters_rf1_multi100_k10_100.RDS")
# filtered ASVs
asvs_rf1_k5 <- readRDS("results/models/asvs_rf1_multi100_k5_10.RDS")
asvs_rf1_k10 <- readRDS("results/models/asvs_rf1_multi100_k10_10.RDS")
# unfiltered_ASVs
asvs_unf_rf2_k5 <- readRDS("results/models/asvs_rf2_multi_k5.RDS")
asvs_unf_rf2_k10 <- readRDS("results/models/asvs_rf2_multi_k10.RDS")


# get k-fold rcv metrics for each set of models
###############################################
# define a function
get_val <-function(modelo, x) {
  y <- lapply(x, function(x) x$trained_model) %>%
    combine_hp_performance() %>%
    `[[`(1) %>%
    mutate(model = modelo)
}

# get results for each model
val_multi <- list()
val_multi[[1]] <- get_val("cl1-k5",cl_rf1_k5)
val_multi[[2]] <- get_val("cl1-k10",cl_rf1_k10) 
val_multi[[3]] <- get_val("asv1-k5",asvs_rf1_k5)
val_multi[[4]] <- get_val("asv1-k10",asvs_rf1_k10)
val_multi[[5]] <- get_val("asvUnF-k5",asvs_unf_rf2_k5)
val_multi[[6]] <- get_val("asvUnF-k10",asvs_unf_rf2_k10)


# make a multi panel plot to compare
p1 <- plot_hp_performance(val_multi[[1]], mtry, AUC) 
p2 <- plot_hp_performance(val_multi[[2]], mtry, AUC) 
p3 <- plot_hp_performance(val_multi[[3]], mtry, AUC) 
p4 <- plot_hp_performance(val_multi[[4]], mtry, AUC) 
p5 <- plot_hp_performance(val_multi[[5]], mtry, AUC) 
p6 <- plot_hp_performance(val_multi[[6]], mtry, AUC) 


ggarrange(p1, p2, p3, p4,p5,p6,
          labels = c("a)", "b)", "c)", "d)", "e)","f)"),
          ncol = 2, nrow = 3)

# get performance metrics for each set of models
###############################################
# define a fucntion
get_performance <-function(modelo, x) {
  lapply(x, function(result) {result[["performance"]] %>% 
  select(cv_metric_AUC, AUC, Accuracy, Kappa, Sensitivity, Specificity)}) %>%
  dplyr::bind_rows() %>%
  mutate(model = modelo)
}

# get performance results for each group of models
perf_multi <- list()

perf_multi[[1]] <- get_performance("cl-k5",cl_rf1_k5)
perf_multi[[2]] <- get_performance("cl-k10",cl_rf1_k10) 
perf_multi[[3]] <- get_performance("asvF-k5",asvs_rf1_k5)
perf_multi[[4]] <- get_performance("asvF-k10",asvs_rf1_k10)
perf_multi[[5]] <- get_performance("asvUnF-k5",asvs_unf_rf2_k5 )
perf_multi[[6]] <- get_performance("asvUnF-k10",asvs_unf_rf2_k10)

# make a boxplot to compare
Reduce(full_join,perf_multi) %>%
  pivot_longer(!model, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=model)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) #+
  facet_wrap(~model)


# Adding SVM model results to compare
#####################################

asvs_svm1_k5 <- readRDS("results/models/asvs_svm1_multi100_k5_10.RDS")

perf_multi[[7]] <- get_performance("asvs-SVM_k5",asvs_svm1_k5)

Reduce(full_join,perf_multi) %>%
  pivot_longer(!model, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x=variable, y=value, fill=model)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) #+
facet_wrap(~model)

