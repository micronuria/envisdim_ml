
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
<<<<<<< HEAD
# 1 - Groups ASVs at Genus level making clusters,
#     reducing biological features from 1,193 to 688.
#     The best taxonomic classification for each cluster
#     is kept even for clusters where genus was not 
#     determined.
# 2 - Add metadata to Clusters
# 3 - Performs clusters data exploration with 
#     proportions (relative abundances) and 
#     rclr transformations.
=======
# 1 - Groups ASVs at Genus level making clusters.
#     The best taxonomic classification for each cluster
#     is kept even for clusters where genus was not 
#     determined.
# 2 - Performs clusters data exploration
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
#---------------------------------------------

############################################
# Load dependencies
############################################

library(tidyverse)
library(vegan)
library(RColorBrewer)
<<<<<<< HEAD
library(viridis)
library(gplots)

library(tictoc)
library(future)
library(doFuture)
library(future.apply)
=======
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893

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

<<<<<<< HEAD
=======
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


>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
# Calculate relative abundances
clusters_ra <- clusters_count %>%
  group_by(sample) %>%                                 
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         
  select(-count)  
  

<<<<<<< HEAD
=======
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

>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
# to wide format
clusters_co_wide <- clusters_count %>%
  pivot_wider(names_from = taxonomy, values_from = count)

clusters_ra_wide <- clusters_ra %>%
  pivot_wider(names_from = taxonomy, values_from = rel_abun)

<<<<<<< HEAD
############################################
# Clusters exploration
############################################

#--------------------------------------
#  Search missing values
#--------------------------------------
table(is.na(clusters_co_wide))
# FALSE 
# 115868  

table(is.na(clusters_ra_wide))
# FALSE 
# 115868 

#--------------------------------------
#  check outcome variable distribution
#--------------------------------------

clusters_co_wide %>%
  with(table(event))
# event
# bloom normal 
# 31    135 

clusters_ra_wide %>%
  with(table(event))
# event
# bloom normal 
# 31    135 

# -------------------------------
# Data distribution
# -------------------------------
# Boxplot count data
pdf(file="results/boxplots/clusters_boxplot_counts.pdf")
clusters_count %>%
  ggplot(aes( y=count,x=taxonomy)) +
  geom_boxplot(outlier.size = 0.4) +
  xlab("Cluster") +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~campaign)
dev.off()

# Boxplot relative abundance data
pdf(file="results/boxplots/clusters_boxplot_ra.pdf")
clusters_ra %>%
  ggplot(aes(y=rel_abun, x=taxonomy)) +
  geom_boxplot(outlier.size = 0.4) +
  xlab("Cluster") +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~campaign)
dev.off()

# Heatmaps
colores <- brewer.pal(n = 2, name = "Dark2")
pdf(file="results/heatmaps/clusters_heatmap_count.pdf")
gplots::heatmap.2(as.matrix(clusters_co_wide[,11:ncol(clusters_co_wide)]),
                  main = "Clusters - count data",
                  xlab = NULL,
                  dendrogram = "none",
                  Colv = FALSE,
                  Rowv = FALSE,
                  trace="none",
                  RowSideColors = colores[clusters_co_wide$campaign],
                  col = viridis(1000))
dev.off()

pdf(file="results/heatmaps/clusters_heatmap_ra.pdf")
gplots::heatmap.2(as.matrix(clusters_ra_wide[,11:ncol(clusters_co_wide)]),
                  main = "Clusters - relative abundance data",
                  xlab = NULL,
                  dendrogram = "none",
                  Colv = FALSE,
                  Rowv = FALSE,
                  trace="none",
                  RowSideColors = colores[clusters_co_wide$campaign],
                  col = viridis(1000))
dev.off()

# ---------------------------------------
# Ordination analysis - NMDS
# ---------------------------------------


=======
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
# define function to plot NMDS results
# coloring samples with different variables
# ---------------------------------------

<<<<<<< HEAD
myPlot <- function(nmds, df_wide) {
  
  
  
  par(mfrow = c(3,3),
      mar = c(4, 4, 2, 1))
=======
myPlot <- function(nmds) {
  
  par(mfrow = c(3,3))
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
  # Stressplot
  stressplot(nmds)
  
  # text
  ordiplot(nmds, display = 'sites', type = 't', main = "sample")
  
  colores = brewer.pal(3, "Set1")
  # Plot - Bloom / no bloom
  ordiplot(nmds, display = 'sites', type = 'n', main ="event")
  points(nmds, pch = 16, cex = 0.9,
<<<<<<< HEAD
         col =colores[as.numeric(as.factor(df_wide$event))])
  legend("topright",legend = as.character(levels(as.factor(df_wide$event))),
         pch = 16, cex = 0.8, col=colores)
=======
         col =colores[as.numeric(as.factor(clusters_co_wide$event))])
  legend("topright",legend = as.character(levels(as.factor(clusters_co_wide$event))),
         pch = 16, cex = 1.0, col=colores)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
  # Plot - Campaign
  ordiplot(nmds, display = 'sites', type = 'n', main="campaign")
  points(nmds, pch = 16, cex = 0.9,
<<<<<<< HEAD
         col =colores[as.numeric(df_wide$campaign)])
  legend("topright",legend = levels(df_wide$campaign),
         pch = 16, cex = 0.8, col=colores)
=======
         col =colores[as.numeric(clusters_co_wide$campaign)])
  legend("topright",legend = levels(clusters_co_wide$campaign),
         pch = 16, cex = 1.0, col=colores)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
  # Plot - site
  ordiplot(nmds, display = 'sites', type = 'n', main = "site")
  points(nmds, pch = 16, cex = 0.9,
<<<<<<< HEAD
         col =colores[df_wide$st])
  legend("topright",legend = as.character(levels(df_wide$st)),
         pch = 16, cex = 0.8, col=colores)
=======
         col =colores[clusters_co_wide$st])
  legend("topright",legend = as.character(levels(clusters_co_wide$st)),
         pch = 16, cex = 1.0, col=colores)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
  
  # Plot - year
  ordiplot(nmds, display = 'sites', type = 'n', main = "year")
  points(nmds, pch = 16, cex = 0.9,
<<<<<<< HEAD
         col =colores[df_wide$year])
  legend("topright",legend = as.character(levels(df_wide$year)),
         pch = 16, cex = 0.8, col=colores)
=======
         col =colores[clusters_co_wide$year])
  legend("topright",legend = as.character(levels(clusters_co_wide$year)),
         pch = 16, cex = 1.0, col=colores)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
  # Plot by month
  colores = colorRampPalette(c("blue","green","orange","red","darkviolet"))(11)
  
<<<<<<< HEAD
  meses <- factor(df_wide$month, 
=======
  meses <- factor(clusters_co_wide$month, 
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
                  levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul','Aug', 'Sep','Oct', 'Nov'))
  
  ordiplot(nmds, display = 'sites', type = 'n', main="month")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[meses])
<<<<<<< HEAD
  legend("topright",legend = as.character(levels(meses)),
         pch = 16, cex = 0.6, col=colores)
=======
  legend("bottomright",legend = as.character(levels(meses)),
         pch = 16, cex = 0.7, col=colores)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
  # Plot by season
  colores = c("orange", "green", "red", "blue")
  ordiplot(nmds, display = 'sites', type = 'n', main="season")
  points(nmds, pch = 16, cex = 0.9,
<<<<<<< HEAD
         col =colores[df_wide$season])
  legend("topright",legend = as.character(levels(df_wide$season)),
         pch = 16, cex = 0.8, col=colores)
=======
         col =colores[clusters_co_wide$season])
  legend("topright",legend = as.character(levels(clusters_co_wide$season)),
         pch = 16, cex = 1, col=colores)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
  # Plot by depth level
  colores = c("blue","green", "red")
  ordiplot(nmds, display = 'sites', type = 'n', main="depth level")
  points(nmds, pch = 16, cex = 0.9,
<<<<<<< HEAD
         col =colores[df_wide$depth_lev])
  legend("topright",legend = levels(df_wide$depth_lev),
         pch = 16, cex = 0.8, col=colores)
  
  par(mfrow = c(1,1),
      mar = c(5, 4, 4, 2)+0.1)
=======
         col =colores[clusters_co_wide$depth_lev])
  legend("topright",legend = levels(clusters_co_wide$depth_lev),
         pch = 16, cex = 1, col=colores)
  
  par(mfrow = c(1,1))
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
  
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
<<<<<<< HEAD
# run in parallel 
nmds.rclr <- metaMDS(euc.d,
                try = 500,
                trymax = 5000,
=======
nmds.rclr <- metaMDS(euc.d,
                try = 500,
                trymax = 1000,
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
                autotransform = FALSE,
                wascores = TRUE,
                expand = FALSE)

<<<<<<< HEAD

# Explore results
# ------------------
nmds.rclr
# Call:
#   metaMDS(comm = euc.d, try = 500, trymax = 5000, autotransform = FALSE,      wascores = TRUE, expand = FALSE) 
=======
# Explore results
# ------------------
nmds.rclr
# 
# Call:
#   metaMDS(comm = euc.d, try = 500, trymax = 1000, autotransform = FALSE,      wascores = TRUE, expand = FALSE) 
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     euc.d 
# Distance: euclidean 
# 
# Dimensions: 2 
<<<<<<< HEAD
# Stress:     0.1404172 
# Stress type 1, weak ties
# Best solution was not repeated after 5000 tries
# The best solution was from try 3708 (random start)
=======
# Stress:     0.1708856 
# Stress type 1, weak ties
# Best solution was repeated 4 times in 500 tries
# The best solution was from try 334 (random start)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
# Scaling: centring, PC rotation 
# Species: scores missing


<<<<<<< HEAD
# NMDS plots
pdf(file = "results/NMDS_plots/clusters_NMDS_counts.pdf")
myPlot(nmds.rclr, clusters_co_wide)
dev.off()
=======

# NMDS plots
myPlot(nmds.rclr)
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893


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
<<<<<<< HEAD
=======

# > nmds.bc
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893
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
<<<<<<< HEAD
# Stress:     0.1390277 
# Stress type 1, weak ties
# Best solution was not repeated after 1000 tries
# The best solution was from try 744 (random start)
# Scaling: centring, PC rotation, halfchange scaling 
# Species: scores missing


# Plots
pdf(file = "results/NMDS_plots/clusters_NMDS_ra.pdf")
myPlot(nmds.bc, clusters_ra_wide)
dev.off()

############################################
# Clean environment
rm(list=c("all_asv", "chlclass","tax_raw","taxonomy", 
          "df_genus.rclr","euc.d", "nmds.rclr","nmds.bc", "df_genus.co",
          "df_genus.ra","tax_raw", "taxonomy","met_raw","myPlot"))
=======
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
>>>>>>> 77630f9a1a2e9e573b2a90e7578e45db40c8d893

# Save objects
save(list=ls(), file="processed_data/clusters_grouping.RData")
