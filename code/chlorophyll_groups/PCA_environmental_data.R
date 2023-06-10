
# Author: Nuria Fernández-Gonzalez

###########################################
#
# PCAs to explore environmental metadata 
#
###########################################

#---------------------------------------------
# Description. This script:
# 1 - PCAs of Envision
# 2 - PCAs of Dimension
# 3 - Barplots for primary production - Dimension
#---------------------------------------------


#######################################################

# PCAs of Envision 16S samples

#######################################################

library(RColorBrewer)

# Metadata preparation
# Load Metadata
metadata <- read.table("./analysis/intermediate/metadata_corrected_16S_table_env.tsv", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
metadata$month <- as.factor(metadata$month)
metadata$st <- as.factor(metadata$st)

rownames(metadata) <- metadata$sample_name

# create a new factor to group depths
source("./analysis/scripts/recode.r")
metadata$depth <- recode(var = metadata$prof_label, "c('P1', 'P2') = 'surface';
                                                   c('P3', 'P4') = 'intermediate'; 
                                                   c('P4', 'P5', 'P6', 'P7') = 'deep'",
                              as.factor = TRUE, levels = c('surface', 'intermediate', 'deep'))


# PCA with variables related to blooms 
mdatabl <- metadata[,c("sample_name","month","day","st","depth","prof",
                       "chlaTotal",
                       "AB","AFSyne","AFProc","AFsp","AFlp",
                       "BB","BFSyne","BFProc","BBsp","BBlp")]

# Remove rows with NAs
mdatabl <- mdatabl[complete.cases(mdatabl),]

# PCA metadata
library(vegan)
# 'species' = columns, 'sites'=rows
PCA <- rda(mdatabl[7:17], scale=TRUE)
summary(PCA)

# variables loadings (correlations to PCs)
loadings <- scores(PCA, display = 'species', scaling = 0)
# correlations of variables to PC1
sort(abs(loadings[,1]), decreasing = TRUE)
# correlations of variables to PC2
sort(abs(loadings[,2]), decreasing = TRUE)

# Plot - Month
colores3 = c("#29a709", "red","blue")
ordiplot(PCA, display = 'sites', type = 'n', 
         main = "Envision - month",
         xlab = "PC1 (48.1 %)",
         ylab = "PC2 (22.3 %)")
points(PCA, pch = 16, cex = 1.3,
       col =colores3[as.numeric(mdatabl$month)])
#orditorp(PCA, display = 'sites')
legend("bottomleft",legend = as.character(levels(mdatabl$month)),
       pch = 16, cex = 1.2, col=colores3)


# Plot - depth
colores3 = c("#0f8135", "#0099ff","#6303bd")
ordiplot(PCA, display = 'sites', type = 'n', 
         main = "Envision - depth",
         xlab = "PC1 (48.1 %)",
         ylab = "PC2 (22.3 %)")
#orditorp(PCA, display = 'sites')
points(PCA, pch = 16, cex = 1.3,
       col =colores3[as.numeric(mdatabl$depth)])
legend("bottomleft",legend = as.character(levels(mdatabl$depth)),
       pch = 16, cex = 1.2, col=colores3)

# Plot - Site
colores2 = c("red","blue")
ordiplot(PCA, display = 'sites', type = 'n',
         main = "Envision - sampling site",
         xlab = "PC1 (48.1 %)",
         ylab = "PC2 (22.3 %)")
points(PCA, pch = 16, cex = 1.3,
       col =colores2[as.numeric(mdatabl$st)])
#orditorp(PCA, display = 'sites')
legend("bottomleft",legend = c("coast-st3","ocean-st6"),
       pch = 16, cex = 1.2, col=colores2)



#######################################################
# PCAs to explore metadata of DIMENSION 16S samples
#######################################################

# Load data
metadim <- read.table("./analysis/intermediate/metadata_corrected_dim.tsv", 
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
metadim$year <- as.factor(metadim$year)
metadim$month <- as.factor(metadim$month)
metadim$Upw <- as.factor(metadim$Upw)
metadim$Sea <- as.factor(metadim$Sea)

# check NAs in numeric variables
lista <- apply(metadim[7:52],2, function(x) {table(is.na(x))})
tablaNAs <- do.call(rbind.data.frame, lista)
names(tablaNAs) <- tolower(names(lista[[5]]))
tablaNAs['true'][tablaNAs['true'] == 36] <- 0
rownames(tablaNAs) <- names(lista)

# Select variables
tokeep <- c("sample","year","month","Dat","Dep","Sea","Upw","Chla.t",
            "PP.p","PP.n","PP.m","BP","PP.t.h")
metadims <- metadim[,tokeep]
rownames(metadims) <- metadim$sample

# Remove rows with NAs
metadims <- metadims[complete.cases(metadims),]

# Calculate PCA
#PCAdim <- rda(metadims[8:15], scale=TRUE)
PCAdim <- rda(metadims[8:13], scale=TRUE)
summary(PCAdim)

# variables loadings (correlations to PCs)
loadingsDIM <- scores(PCAdim, display = 'species', scaling = 0)
# correlations of variables to PC1
sort(abs(loadingsDIM[,1]), decreasing = TRUE)
# correlations of variables to PC2
sort(abs(loadingsDIM[,2]), decreasing = TRUE)

# Plot - month
colores = c(brewer.pal(9, "Set1"), brewer.pal(2, "Set2"))
ordiplot(PCAdim, display = 'sites', type = 'n', 
         main = "Dimension - month",
         xlab = "PC1 (63.1 %)",
         ylab = "PC2 (18.2 %)")
points(PCAdim, pch = 16, cex = 1.4,
       col =colores[as.numeric(metadims$month)])
legend("bottomright",legend = as.character(levels(metadims$month)),
       pch = 16, cex = 1, col=colores)

# Plot - season
colores4 = c("blue","green","red","#f78215")
ordiplot(PCAdim, display = 'sites', type = 'n',
         main = "Dimension - season",
         xlab = "PC1 (63.1 %)",
         ylab = "PC2 (18.2 %)")
points(PCAdim, pch = 16, cex = 1.3,
       col =colores4[as.numeric(metadims$Sea)])
legend("bottomright",legend = as.character(levels(metadims$Sea)),
       pch = 16, cex = 1.2, col=colores4)

# Plot - depth
colores2 = c("red","blue")
ordiplot(PCAdim, display = 'sites', type = 'n',
         main = "Dimension - depth (m)",
         xlab = "PC1 (63.1 %)",
         ylab = "PC2 (18.2 %)")
points(PCAdim, pch = 16, cex = 1.3,
       col =colores2[as.numeric(as.factor(metadims$Dep))])
legend("bottomright",legend = as.character(unique(metadims$Dep)),
       pch = 16, cex = 1.2, col=colores2)

# Plot - year
colores2 = c("red","blue")
ordiplot(PCAdim, display = 'sites', type = 'n',
         main = "Dimension - year",
         xlab = "PC1 (63.1 %)",
         ylab = "PC2 (18.2 %)")
points(PCAdim, pch = 16, cex = 1.3,
       col =colores2[as.numeric(as.factor(metadims$year))])
legend("bottomright",legend = as.character(unique(metadims$year)),
       pch = 16, cex = 1.2, col=colores2)


# Plot - upwelling
colores2 = c("red","blue")
ordiplot(PCAdim, display = 'sites', type = 'n',
         main = "Dimension - upwelling events",
         xlab = "PC1 (63.1 %)",
         ylab = "PC2 (18.2 %)")
points(PCAdim, pch = 16, cex = 1.3,
       col =colores2[as.numeric(as.factor(metadims$Upw))])
legend("bottomright",legend = as.character(unique(metadims$Upw)),
       pch = 16, cex = 1.2, col=colores2)

#######################################################
# Barplot of primary production of Dimension dataset
#######################################################

# Create a new variable to plot all samples
metadims$dia <- paste(metadims$Dat, metadims$Dep, sep = "-")

euk <- as_labeller(
  c(`PP.n` = "Nanoeukaryotes",
    `PP.m` = "Microeukaryotes",
    `PP.p` = "Picoeukaryotes"))

# barplots
tibble(metadims) %>% 
  select(-PP.t.h) %>%
  pivot_longer(cols = starts_with("PP"),
               names_to = "type", values_to = "PrimP") %>%
  ggplot(., aes(x=fct_inorder(dia), y=PrimP, fill = month)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = colores) +
    labs(x="Date - Depth", y = "Primary production (mgC/m3h)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_wrap(vars(type), nrow = 3, labeller = euk)

