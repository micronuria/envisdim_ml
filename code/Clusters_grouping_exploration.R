
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# DATA PREPROCESS
# grouping ASVs to clusters
#
###########################################
#
#---------------------------------------------
# Description. This script:
# 1 - Groups ASVs at Genus level making clusters,
#     reducing biological features from 1,193 to 688.
#     The best taxonomic classification for each cluster
#     is kept even for clusters where genus was not 
#     determined.
# 2 - Add metadata to Clusters
# 3 - Performs clusters data exploration with 
#     proportions (relative abundances) and 
#     rclr transformations.
#---------------------------------------------

############################################
# Load dependencies
############################################

library(tidyverse)
library(vegan)
library(RColorBrewer)
library(viridis)
library(gplots)

library(tictoc)
library(future)
library(doFuture)
library(future.apply)

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

# Calculate relative abundances
clusters_ra <- clusters_count %>%
  group_by(sample) %>%                                 
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         
  select(-count)  
  

# to wide format
clusters_co_wide <- clusters_count %>%
  pivot_wider(names_from = taxonomy, values_from = count)

clusters_ra_wide <- clusters_ra %>%
  pivot_wider(names_from = taxonomy, values_from = rel_abun)

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


# define function to plot NMDS results
# coloring samples with different variables
# ---------------------------------------

myPlot <- function(nmds, df_wide) {
  
  
  
  par(mfrow = c(3,3),
      mar = c(4, 4, 2, 1))
  
  # Stressplot
  stressplot(nmds)
  
  # text
  ordiplot(nmds, display = 'sites', type = 't', main = "sample")
  
  colores = brewer.pal(3, "Set1")
  # Plot - Bloom / no bloom
  ordiplot(nmds, display = 'sites', type = 'n', main ="event")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[as.numeric(as.factor(df_wide$event))])
  legend("topright",legend = as.character(levels(as.factor(df_wide$event))),
         pch = 16, cex = 0.8, col=colores)
  
  # Plot - Campaign
  ordiplot(nmds, display = 'sites', type = 'n', main="campaign")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[as.numeric(df_wide$campaign)])
  legend("topright",legend = levels(df_wide$campaign),
         pch = 16, cex = 0.8, col=colores)
  
  # Plot - site
  ordiplot(nmds, display = 'sites', type = 'n', main = "site")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[df_wide$st])
  legend("topright",legend = as.character(levels(df_wide$st)),
         pch = 16, cex = 0.8, col=colores)
  
  
  # Plot - year
  ordiplot(nmds, display = 'sites', type = 'n', main = "year")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[df_wide$year])
  legend("topright",legend = as.character(levels(df_wide$year)),
         pch = 16, cex = 0.8, col=colores)
  
  # Plot by month
  colores = colorRampPalette(c("blue","green","orange","red","darkviolet"))(11)
  
  meses <- factor(df_wide$month, 
                  levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul','Aug', 'Sep','Oct', 'Nov'))
  
  ordiplot(nmds, display = 'sites', type = 'n', main="month")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[meses])
  legend("topright",legend = as.character(levels(meses)),
         pch = 16, cex = 0.6, col=colores)
  
  # Plot by season
  colores = c("orange", "green", "red", "blue")
  ordiplot(nmds, display = 'sites', type = 'n', main="season")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[df_wide$season])
  legend("topright",legend = as.character(levels(df_wide$season)),
         pch = 16, cex = 0.8, col=colores)
  
  # Plot by depth level
  colores = c("blue","green", "red")
  ordiplot(nmds, display = 'sites', type = 'n', main="depth level")
  points(nmds, pch = 16, cex = 0.9,
         col =colores[df_wide$depth_lev])
  legend("topright",legend = levels(df_wide$depth_lev),
         pch = 16, cex = 0.8, col=colores)
  
  par(mfrow = c(1,1),
      mar = c(5, 4, 4, 2)+0.1)
  
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
# run in parallel 
nmds.rclr <- metaMDS(euc.d,
                try = 500,
                trymax = 5000,
                autotransform = FALSE,
                wascores = TRUE,
                expand = FALSE)


# Explore results
# ------------------
nmds.rclr
# Call:
#   metaMDS(comm = euc.d, try = 500, trymax = 5000, autotransform = FALSE,      wascores = TRUE, expand = FALSE) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     euc.d 
# Distance: euclidean 
# 
# Dimensions: 2 
# Stress:     0.1404172 
# Stress type 1, weak ties
# Best solution was not repeated after 5000 tries
# The best solution was from try 3708 (random start)
# Scaling: centring, PC rotation 
# Species: scores missing


# NMDS plots
pdf(file = "results/NMDS_plots/clusters_NMDS_counts.pdf")
myPlot(nmds.rclr, clusters_co_wide)
dev.off()


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

# Save objects
save(list=ls(), file="processed_data/clusters_grouping.RData")
