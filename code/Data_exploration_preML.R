# Load preprocessed data 
source('code/ASVs_preprocess.R')

library(vegan)

# Prepare data to work with vegan package 
df_asvs <- as.data.frame(composite_wide)
rownames(df_asvs)<-df_asvs$sample
evento <- df_asvs$event
df_asvs <- df_asvs[,4:ncol(df_asvs)]

# apply robust center log ratio transformation to count data
df_asvs.rclr <- decostand(df_asvs, "rclr", MARGIN = 1)

# Calculate euclidean distances
euc.d <- vegdist(df_asvs.rclr, method="euclidean")

# Calculate a non-metric multidimensional analysis
nmds <- metaMDS(euc.d,
                try = 500,
                trymax = 500,
                autotransform = FALSE,
                wascores = TRUE,
                expand = FALSE)
nmds

# Screeplot
stressplot (nmds)

# Plot
colores = c("red","blue")
ordiplot(nmds, display = 'sites', type = 't')
points(nmds, pch = 16, cex = 0.7,
       col =colores[as.numeric(as.factor(evento))])
legend("topright",legend = as.character(levels(as.factor(evento))),
       pch = 16, cex = 1.0, col=colores)



