
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# DATA PREPROCESS
# grouping ASVs to clusters
#
###########################################

#---------------------------------------------
# Description. This script:
# 1 - Groups ASVs at Genus level making clusters.
#     The best taxonomic classification for each cluster
#     is kept even for clusters where genus was not 
#     determined.
# 2 - Performs clusters data exploration
#---------------------------------------------

############################################
# Load dependencies
############################################

library(tidyverse)
library(vegan)
library(RColorBrewer)

set.seed(20210401)

load("processed_data/initial_preprocess.RData")

############################################
# Group data to Genus level
############################################

# Group taxonomy to filter data and format taxonomy
# to get the last taxonomic level assignment.
taxonomy <- tax_raw %>%
  select(asv, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  unite("taxonomy", Kingdom:Genus, sep = ";", remove = TRUE) %>%
  mutate(taxonomy = str_replace(taxonomy, ";NA","_unclassified"),
         taxonomy = str_replace_all(taxonomy, ";NA", ""),
         taxonomy = str_replace_all(taxonomy, ".*;", ""))

# join count data, taxonomy, metadata, chlorophyll groups 
# and change some variables to factors
clusters_count <- inner_join(all_asv, taxonomy, by = "asv") %>%
  group_by(sample, taxonomy) %>%                        # group by genus
  summarize(count = sum(count), .groups="drop") %>%     # get rid of count column
  inner_join(., chlclass, by="sample") %>%              # add data about chlorophyll groups
  inner_join(., met_raw) %>%
  mutate(gr_event = as.factor(gr_event),
         campaign = as.factor(campaign),
         st = as.factor(st),
         year = as.factor(year),
         month = as.factor(month),
         season = as.factor(season),
         depth_lev = as.factor(depth_lev))

# > clusters_count
# # A tibble: 114,208 × 12
# sample   taxonomy                             count event  campaign  gr_event   st    year  month depth season depth_lev
# <chr>    <chr>                                <int> <chr>  <fct>     <fct>      <fct> <fct> <fct> <dbl> <fct>  <fct>    
#   1 DIMProk1 028H05-P-BN-P5_unclassified              0 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 2 DIMProk1 053A03-B-DI-P58_unclassified            57 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 3 DIMProk1 211ds20_unclassified                     0 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 4 DIMProk1 A714019                                571 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 5 DIMProk1 AB1_unclassified                         0 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 6 DIMProk1 AEGEAN-169 marine group_unclassified  1516 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 7 DIMProk1 AT-s2-59_unclassified                    0 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 8 DIMProk1 AT-s3-44_unclassified                  257 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 9 DIMProk1 Acanthopleuribacter                     28 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# 10 DIMProk1 Acetobacteraceae_unclassified            2 normal dimension normal_dim 3     2014  Jan       5 Winter surface  
# # ℹ 114,198 more rows
# # ℹ Use `print(n = ...)` to see more rows


# Calculate relative abundances
clusters_ra <- clusters_count %>%
  group_by(sample) %>%                                 
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         
  select(-count)  
  

############################################
# Clusters exploration
############################################

# -------------------------------
# Data distribution
# -------------------------------
# count data
clusters_count %>%
  ggplot(aes(count, taxonomy)) +
  geom_boxplot() +
  ylab("Cluster") +
  theme(axis.text.y = element_blank())

# relative abundance data
clusters_ra %>%
  ggplot(aes(rel_abun, taxonomy)) +
  geom_boxplot() +
  ylab("Cluster") +
  theme(axis.text.y = element_blank())

# ---------------------------------------
# Ordination analysis - NMDS
# ---------------------------------------

# to wide format
clusters_co_wide <- clusters_count %>%
  pivot_wider(names_from = taxonomy, values_from = count)

clusters_ra_wide <- clusters_ra %>%
  pivot_wider(names_from = taxonomy, values_from = rel_abun)

# define function to plot NMDS results
# coloring samples with different variables
# ---------------------------------------

myPlot <- function(nmds) {
  
  par(mfrow = c(3,3))
  
  # Stressplot
  stressplot(nmds)
  
  # text
  ordiplot(nmds, display = 'sites', type = 't', main = "sample")
  
  colores = brewer.pal(3, "Set1")
  # Plot - Bloom / no bloom
  ordiplot(nmds, display = 'sites', type = 'n', main ="event")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[as.numeric(as.factor(clusters_co_wide$event))])
  legend("topright",legend = as.character(levels(as.factor(clusters_co_wide$event))),
         pch = 16, cex = 1.0, col=colores)
  
  # Plot - Campaign
  ordiplot(nmds, display = 'sites', type = 'n', main="campaign")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[as.numeric(clusters_co_wide$campaign)])
  legend("topright",legend = levels(clusters_co_wide$campaign),
         pch = 16, cex = 1.0, col=colores)
  
  # Plot - site
  ordiplot(nmds, display = 'sites', type = 'n', main = "site")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[clusters_co_wide$st])
  legend("topright",legend = as.character(levels(clusters_co_wide$st)),
         pch = 16, cex = 1.0, col=colores)
  
  
  # Plot - year
  ordiplot(nmds, display = 'sites', type = 'n', main = "year")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[clusters_co_wide$year])
  legend("topright",legend = as.character(levels(clusters_co_wide$year)),
         pch = 16, cex = 1.0, col=colores)
  
  # Plot by month
  colores = colorRampPalette(c("blue","green","orange","red","darkviolet"))(11)
  
  meses <- factor(clusters_co_wide$month, 
                  levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul','Aug', 'Sep','Oct', 'Nov'))
  
  ordiplot(nmds, display = 'sites', type = 'n', main="month")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[meses])
  legend("bottomright",legend = as.character(levels(meses)),
         pch = 16, cex = 0.7, col=colores)
  
  # Plot by season
  colores = c("orange", "green", "red", "blue")
  ordiplot(nmds, display = 'sites', type = 'n', main="season")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[clusters_co_wide$season])
  legend("topright",legend = as.character(levels(clusters_co_wide$season)),
         pch = 16, cex = 1, col=colores)
  
  # Plot by depth level
  colores = c("blue","green", "red")
  ordiplot(nmds, display = 'sites', type = 'n', main="depth level")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[clusters_co_wide$depth_lev])
  legend("topright",legend = levels(clusters_co_wide$depth_lev),
         pch = 16, cex = 1, col=colores)
  
  par(mfrow = c(1,1))
  
}


# ---------------------------------------
# calculate NMDS - rclr transformation
# ---------------------------------------

# prepare data
df_genus.co <- clusters_co_wide[,11:ncol(clusters_co_wide)]
rownames(df_genus.co)<-clusters_co_wide$sample

# apply robust center log ratio transformation to count data
df_genus.rclr <- decostand(df_genus.co, "rclr", MARGIN = 1)

# Calculate euclidean distances
euc.d <- vegdist(df_genus.rclr, method="euclidean")

# Calculate a non-metric multidimensional analysis
nmds.rclr <- metaMDS(euc.d,
                try = 500,
                trymax = 1000,
                autotransform = FALSE,
                wascores = TRUE,
                expand = FALSE)

# Explore results
# ------------------
nmds.rclr
# 
# Call:
#   metaMDS(comm = euc.d, try = 500, trymax = 1000, autotransform = FALSE,      wascores = TRUE, expand = FALSE) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     euc.d 
# Distance: euclidean 
# 
# Dimensions: 2 
# Stress:     0.1708856 
# Stress type 1, weak ties
# Best solution was repeated 4 times in 500 tries
# The best solution was from try 334 (random start)
# Scaling: centring, PC rotation 
# Species: scores missing



# NMDS plots
myPlot(nmds.rclr)


# ---------------------------------------------
# calculate NMDS - with Bray-Curtis dissimilarities
# ---------------------------------------------

# prepare data
df_genus.ra <- clusters_ra_wide[,11:ncol(clusters_ra_wide)]
rownames(df_genus.ra)<-clusters_ra_wide$sample


# Calculate a non-metric multidimensional analysis
nmds.bc <- metaMDS(df_genus.ra,
                   distance = "bray",
                   try = 500,
                   trymax = 1000,
                   autotransform = TRUE,
                   wascores = FALSE,
                   expand = FALSE)

# Explore results
# ------------------
nmds.bc

# > nmds.bc
# 
# Call:
#   metaMDS(comm = df_genus.ra, distance = "bray", try = 500, trymax = 1000,      autotransform = TRUE, wascores = FALSE, expand = FALSE) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     df_genus.ra 
# Distance: bray 
# 
# Dimensions: 2 
# Stress:     0.1694258 
# Stress type 1, weak ties
# Best solution was repeated 2 times in 500 tries
# The best solution was from try 274 (random start)
# Scaling: centring, PC rotation, halfchange scaling 
# Species: scores missing

# Plots
myPlot(nmds.bc)


############################################
# Clean environment
rm(list=c("all_asv", "chlclass","tax_raw","taxonomy", "df_genus", 
          "df_genus.rclr","euc.d", "nmds.rclr","nmds.bc", "df_genus.co",
          "df_genus.ra", "df_genus.rclr","tax_raw", "taxonomy","met_raw"))

# Save objects
save(list=ls(), file="processed_data/clusters_grouping.RData")
