
# Author: Nuria Fernández-Gonzalez

###########################################
#
# Get feature importance
#
###########################################

#---------------------------------------------
# Description. This script obtains the
# scaled feature importances of the RF model
# trained with filtered ASVs implemented with 
# method SMOTE using the
# varImp function from caret
#---------------------------------------------

############################################
# Load dependencies
############################################

# Load libraries
library(tidyverse)
library(mikropml)

# load models 
results <- readRDS("results/models/test_syntetic_results.RDS")

###############################################
# Define function to extract feature importance
###############################################

get_importance <- function(modelo) {
  
  model <- modelo %>%
    pluck("trained_model")
  
  caret::varImp(model, scale = TRUE, useModel = TRUE) %>% 
    pluck("importance") %>%
    as_tibble(rownames = "feature") %>%
    mutate(value = Overall, .keep = "unused") 
}

# Iterate over all models
var_imp_results <- map_dfr(results, get_importance)


# Make a plot
var_imp_results %>%
  group_by(feature) %>%
  summarize(median = median(value),
            l_quartile = quantile(value, probs = 0.25),
            u_quartile = quantile(value, probs = 0.75)) %>%
  mutate(feature = fct_reorder(feature, median)) %>%
  filter(median > 10) %>%
  
  ggplot(aes(x = median, y = feature, xmin=l_quartile, xmax=u_quartile)) +
    geom_point()+
    geom_linerange()
     
     
ggsave("results/features/importance_over10.tiff", width = 5, height = 8)





