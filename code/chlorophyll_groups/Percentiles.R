
# Author: Nuria Fernández-Gonzalez

#######################################################################
#
#  Total Chlorophyl Quartiles Data Exploration
#  
#######################################################################


#---------------------------------------------
# Description. This script:
# 1 - Data divisions with 50 and 75 percentiles
# 2 - Data divisions with 90 percentile for
#     ocean station
#---------------------------------------------



# Division of data based on 50 and 75 percentiles.
#######################################################################
# Divide chlorophyl data into groups based on different percentils
# to find which samples belong to blooms and which do not.

# Load metadata of Envision campaign
rawenv <- read.table("../intermediate/metadata_corrected_16S_table_env.tsv", 
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(rawenv) <- rawenv$sample_name
# Load metadata of Dimension campaign
rawdim <- read.table("../intermediate/metadata_corrected_dim.tsv",
                     sep ="\t", header = TRUE, stringsAsFactors = FALSE)
rownames(rawdim) <- rawdim$sample

# Separate Envision in the two sations
envst3 <- subset(rawenv, st == 3)
envst6 <- subset(rawenv, st == 6)

# Vectors of chlorophyl value for oceanic and coastal stations
chl_oce <- envst6$chlaTotal
names(chl_oce) <- rownames(envst6)
chl_cos <- c(envst3$chlaTotal, rawdim$Chla.t)
names(chl_cos) <- c(rownames(envst3), rownames(rawdim))

# get statistics
sum_oce <- summary(chl_oce)
sum_cos <- summary(chl_cos)

# Divisions based on the 50 percentile 
bloom_oce_50 <- chl_oce[chl_oce >= sum_oce[3]]
nobloom_oce_50 <- chl_oce[chl_oce < sum_oce[3]]

bloom_cos_50 <- chl_cos[chl_cos >= sum_cos[3]]
nobloom_cos_50 <- chl_cos[chl_cos < sum_cos[3]]

# Divisions based on the 75 percentile
bloom_oce_75 <- chl_oce[chl_oce >= sum_oce[5]]
nobloom_oce_75 <- chl_oce[chl_oce < sum_oce[5]]

bloom_cos_75 <- chl_cos[chl_cos >= sum_cos[5]]
nobloom_cos_75 <- chl_cos[chl_cos < sum_cos[5]]

valores <- list(bloom_cos_50, bloom_cos_75, nobloom_cos_50, nobloom_cos_75,
                bloom_oce_50, bloom_oce_75, nobloom_oce_50, nobloom_oce_75)

names(valores)<- c("bl_C_50","bl_C_75","Nbl.C.50","Nbl_C_75",
                   "bl_O_50","bl_O_75", "Nbl_O_50", "Nbl_O_75")

# Check groups
barplot(unlist(lapply(valores,length)), main = "Number of samples on each group",
        ylab = "frequency", las = 2)

# Values of total Chrlorophyl on each group
# boxplot
boxplot(valores, ylab = "Total Chlorophyl", las = 2)

## Details on the 75 percentile groups
# Check the sample distributions through months and depths for this division.

# Ocean station - ENVISION

df_bloom_oce_75 <- envst6[envst6$sample_name %in% names(bloom_oce_75),]
df_nobloom_oce_75 <- envst6[envst6$sample_name %in% names(nobloom_oce_75),]

# Month
par(mfrow=c(1,2))
barplot(table(df_bloom_oce_75$month), ylab = "frequency", main="bloom", ylim = c(0, 20))
barplot(table(df_nobloom_oce_75$month), ylab = "frequency", main="no bloom", ylim = c(0, 20))

# Depth
par(mfrow=c(1,2))
hist(df_bloom_oce_75$prof, xlab = "Depth (m)", ylab = "frequency", main = "bloom", ylim = c(0,12))
hist(df_nobloom_oce_75$prof,  xlab = "Depth (m)", ylab = "frequency", main = "no bloom",  ylim = c(0,12))

#Coast - ENVISION
df_bloom_cosenv_75 <- envst3[envst3$sample_name %in% names(bloom_cos_75),]

# Month
par(mfrow=c(1,2))
barplot(table(df_bloom_cosenv_75$month), ylab = "frequency", main="bloom", ylim = c(0, 20))
barplot(table(df_nobloom_cosenv_75$month), ylab = "frequency", main="no bloom", ylim = c(0, 20))

# Depth*
par(mfrow=c(1,2))
hist(df_bloom_cosenv_75$prof, xlab = "Depth (m)", ylab = "frequency", main = "bloom", ylim = c(0,12))
hist(df_nobloom_cosenv_75$prof,  xlab = "Depth (m)", ylab = "frequency", main = "no bloom",  ylim = c(0,12))

# Coast - DIMENSION
df_bloom_cosdim_75 <- rawdim[rawdim$sample %in% names(bloom_cos_75),]
df_nobloom_cosdim_75 <- rawdim[rawdim$sample %in% names(nobloom_cos_75),]

# Month
par(mfrow=c(1,2))
barplot(table(df_bloom_cosdim_75$month), ylab = "frequency", main="bloom", ylim = c(0, 5), las =2)
barplot(table(df_nobloom_cosdim_75$month), ylab = "frequency", main="no bloom", ylim = c(0, 5), las =2 )

# Season
par(mfrow=c(1,2))
barplot(table(df_bloom_cosdim_75$Sea), ylab = "frequency", main="bloom", ylim = c(0, 10), las =2)
barplot(table(df_nobloom_cosdim_75$Sea), ylab = "frequency", main="no bloom", ylim = c(0, 10), las =2 )

# Depth
par(mfrow=c(1,2))
hist(df_bloom_cosdim_75$Dep, xlab = "Depth (m)", ylab = "frequency", main = "bloom", ylim = c(0,16))
hist(df_nobloom_cosdim_75$Dep,  xlab = "Depth (m)", ylab = "frequency", main = "no bloom",  ylim = c(0,16))


## 90 percentile for oceanic station
# get the 90 percentile of chlorophyl concentrations on the oceanic station
oce90 <-quantile(chl_oce, probs = c(0.90))

# Divisions based on the 90 percentile 
bloom_oce_90 <- chl_oce[chl_oce >= oce90]
nobloom_oce_90 <- chl_oce[chl_oce < oce90]

valores <- list(bloom_oce_90, nobloom_oce_90)
names(valores)<- c("bl_C_90", "Nbl_C_90")

# check groups
barplot(unlist(lapply(valores,length)), main = "Number of samples on each group",
        ylab = "frequency", las = 2)

# data selection
df_bloom_oce_90 <- envst6[envst6$sample_name %in% names(bloom_oce_90),]
df_nobloom_oce_90 <- envst6[envst6$sample_name %in% names(nobloom_oce_90),]

# Month
par(mfrow=c(1,2))
barplot(table(df_bloom_oce_90$month), ylab = "frequency", main="bloom", ylim = c(0, 25))
barplot(table(df_nobloom_oce_90$month), ylab = "frequency", main="no bloom", ylim = c(0, 25))

# Depth*
par(mfrow=c(1,2))
hist(df_bloom_oce_90$prof, xlab = "Depth (m)", ylab = "frequency", main = "bloom", ylim = c(0,15))
hist(df_nobloom_oce_90$prof,  xlab = "Depth (m)", ylab = "frequency", main = "no bloom",  ylim = c(0,15))

