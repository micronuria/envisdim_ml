
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# DATA PREPROCESS
# ASVs to clusters
#
###########################################

#---------------------------------------------
# Description. This script:
# 1 - Filters ASVs with mean relative abundances < 0.01%
#     that reduces ASVs from 7,569 to 1,193
# 2 - Formats taxonomy and adds metadata
# 3 - Performs ASVs data exploration
#---------------------------------------------

############################################
# Load dependencies
############################################

library(tidyverse)
library(vegan)
library(RColorBrewer)

set.seed(20210401)

load("processed_data/initial_preprocess.RData")


##############################
# pre-filter ASVs
##############################

# Revome low abundance ASVS (below 0.01%)
filtered_asv <- all_asv %>%
  group_by(sample) %>%                                 
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         
  select(-count) %>%
  pivot_wider(names_from = sample, values_from = rel_abun) %>%
  rowwise %>%
  mutate(mean = mean(c_across(c(DIMProk1,ENVProk178))))%>%
  filter(mean > 0.0001) %>% 
  select(-mean) %>%
  pivot_longer(!asv, names_to = "sample", values_to="rel_abun")


############################################
# Format taxonomic assignment to ASVs
############################################

#  format taxonomy to get the last taxonomic level assignment.
tax_asvs <- tax_raw %>%
  filter(asv %in% unique(filtered_asv$asv)) %>%
   unite("taxonomy", Kingdom:Species, sep = ";", remove = TRUE) %>%
   mutate(taxonomy = str_replace_all(taxonomy, "NA","unclassified"))

# join asvs count data, metadata, chlorophyll groups 
# and change some variables to factors.
# Not adding taxonomy as it produces problems with pivot_wider

# relative abuncances
asvs_ra <- inner_join(filtered_asv, chlclass, by="sample")  %>%              
   inner_join(., met_raw, by="sample") %>%
   mutate(gr_event = as.factor(gr_event),
         campaign = as.factor(campaign),
         st = as.factor(st),
         year = as.factor(year),
         month = as.factor(month),
         season = as.factor(season),
         depth_lev = as.factor(depth_lev))

# with count data
asvs_co <- all_asv %>%
  filter(asv %in% unique(filtered_asv$asv)) %>%
  inner_join(., chlclass, by="sample")  %>%
  inner_join(., met_raw, by="sample") %>%
  mutate(gr_event = as.factor(gr_event),
         campaign = as.factor(campaign),
         st = as.factor(st),
         year = as.factor(year),
         month = as.factor(month),
         season = as.factor(season),
         depth_lev = as.factor(depth_lev))


############################################
# ASVs exploration
############################################

# -------------------------------
# Data distribution
# -------------------------------
# count data
pdf(file="results/boxplots/asvs_boxplot_counts.pdf")
asvs_co %>%
  ggplot(aes( y=count,x=asv)) +
  geom_boxplot(outlier.size = 0.4) +
  xlab("ASV") +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~campaign)
dev.off()


# relative abundance data
pdf(file="results/boxplots/asvs_boxplots_relab.pdf")
asvs_ra %>%
  ggplot(aes(y=rel_abun, x=asv)) +
  geom_boxplot(outlier.size = 0.4) +
  xlab("ASV") +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~campaign)
dev.off()


# ---------------------------------------
# Ordination analysis - NMDS
# ---------------------------------------

# to wide format
asvs_co_wide <- asvs_co %>%
  pivot_wider(names_from = asv, values_from = count)

asvs_ra_wide <- asvs_ra %>% 
  pivot_wider(names_from = asv, values_from = rel_abun)

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
df_asvs.co <- asvs_co_wide[,11:ncol(asvs_co_wide)]
rownames(df_asvs.co)<-asvs_co_wide$sample

# apply robust center log ratio transformation to count data
df_asvs.rclr <- decostand(df_asvs.co, "rclr", MARGIN = 1)

# Calculate euclidean distances
euc.d <- vegdist(df_asvs.rclr, method="euclidean")

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

# Call:
#   metaMDS(comm = euc.d, try = 500, trymax = 1000, autotransform = FALSE,      wascores = TRUE, expand = FALSE) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     euc.d 
# Distance: euclidean 
# 
# Dimensions: 2 
# Stress:     0.1579144 
# Stress type 1, weak ties
# Best solution was repeated 3 times in 500 tries
# The best solution was from try 105 (random start)
# Scaling: centring, PC rotation 
# Species: scores missing

# NMDS plots
pdf(file="results/NMDS_plots/asvs_NMDS_counts.pdf")
myPlot(nmds.rclr, asvs_co_wide)
dev.off()

# ---------------------------------------------
# calculate NMDS - with Bray-Curtis dissimilarities
# ---------------------------------------------

# prepare data
df_asvs.ra <- asvs_ra_wide[,11:ncol(asvs_ra_wide)]
rownames(df_asvs.ra)<-asvs_ra_wide$sample


# Calculate a non-metric multidimensional analysis
nmds.bc <- metaMDS(df_asvs.ra,
                   distance = "bray",
                   try = 500,
                   trymax = 1000,
                   autotransform = TRUE,
                   wascores = FALSE,
                   expand = FALSE)

# Explore results
nmds.bc
# Call:
#   metaMDS(comm = df_asvs.ra, distance = "bray", try = 500, trymax = 1000,      autotransform = TRUE, wascores = FALSE, expand = FALSE) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     df_asvs.ra 
# Distance: bray 
# 
# Dimensions: 2 
# Stress:     0.1409043 
# Stress type 1, weak ties
# Best solution was repeated 2 times in 500 tries
# The best solution was from try 371 (random start)
# Scaling: centring, PC rotation, halfchange scaling 
# Species: scores missing

# plot
pdf(file="results/NMDS_plots/asvs_NMDS_ra.pdf")
myPlot(nmds.bc, asvs_ra_wide)
dev.off()


###############################
# clean environment
rm("filtered_asv","met_raw","chlclass","all_asv","tax_raw",
   "df_asvs.rclr","euc.d", "nmds.rclr","df_asvs.co", 
   "df_asvs.ra", "nmds.bc","myPlot")

# Save objects
save(list=ls(), file="processed_data/asvs_filtered.RData")
