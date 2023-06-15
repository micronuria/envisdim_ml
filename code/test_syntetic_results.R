library(tidyverse)

results <- readRDS("results/models/test_syntetic_results.RDS")


prueba <- results[[24]]$performance

# sacar valores de k-fold CV para comparar overfitting

get_performance <-function(x) {
  
  lapply(x, function(result) {result[["performance"]] %>% 
      select(cv_metric_AUC, AUC, Accuracy, Kappa, Sensitivity, Specificity)}) %>%
    dplyr::bind_rows() %>%
    mutate(model = modelo)
}


get_perf <- function(x) {
  a <- lapply(x, function(result){ result$performance$overall[c(1,2)]})
  b <- lapply(x, function(result){ result$performance$byClass[c(1,2)]})
  
}


dfa <- data.frame(matrix(unlist(a), nrow=length(a), byrow=TRUE))
dfb <- data.frame(matrix(unlist(b), nrow=length(b), byrow=TRUE))
prueba <- cbind(dfa,dfb)
names(prueba) <- c(names(a[[1]]), names(b[[1]]))





