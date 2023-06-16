
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

# load taxonomy table
taxa <- read.table("raw_data/asvs_taxonomy.tsv", sep = "\t")
names(taxa) <- taxa[1,]
taxa <- taxa[-1,]

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
  filter(median > 20) %>%
    ggplot(aes(x = median, y = feature, xmin=l_quartile, xmax=u_quartile)) +
    geom_point()+
    geom_linerange()
     
# save plot
ggsave("results/features/importance_over20.tiff", width = 5, height = 8)


# get taxonomic assigments of ASVs
asvs_names <- var_imp_results %>%
  group_by(feature) %>%
  summarize(median = median(value),
            l_quartile = quantile(value, probs = 0.25),
            u_quartile = quantile(value, probs = 0.75)) %>%
  mutate(feature = fct_reorder(feature, median)) %>%
  filter(median > 20) %>%
  select(feature) %>%
  filter(str_detect(feature, "^asv"))

asvs_tax <- subset(taxa, taxa$asv %in% asvs_names$feature)
write.csv2(asvs_tax, file = "results/features/asvs_important_taxa.csv")
