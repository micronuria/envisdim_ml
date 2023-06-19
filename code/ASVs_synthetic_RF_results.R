
# Author: Nuria Fernández-Gonzalez

###########################################
#
# Get results of RF model with synthetic data
#
###########################################

#---------------------------------------------
# Description. This script gets the performance
# of the RF models using synthetic data
#---------------------------------------------

############################################
# Load dependencies
############################################

library(tidyverse)

# load the 100 models
results <- readRDS("results/models/test_syntetic_results.RDS")

# get values for Accuracy and Kappa
a <- lapply(results, function(result){ result$performance$overall[c(1,2)]})
# get values for Sensitivity and Specificity
b <- lapply(results, function(result){ result$performance$byClass[c(1,2)]})

# join in a single dataframe
dfa <- data.frame(matrix(unlist(a), nrow=length(a), byrow=TRUE))
dfb <- data.frame(matrix(unlist(b), nrow=length(b), byrow=TRUE))
performance <- cbind(dfa,dfb)
names(performance) <- c(names(a[[1]]), names(b[[1]]))


# plot
performance %>% as_tibble() %>%
  pivot_longer(everything(), values_to = "value", names_to = "variable") %>%
  ggplot(aes(x=variable, y=value)) + 
  geom_boxplot() + 
  ylim(0, 1) + 
  theme(axis.title = element_text(size = 14))


