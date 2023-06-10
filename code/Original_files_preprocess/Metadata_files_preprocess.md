
# Environmental metadata files preprocess

Data gathered from collaborators lack a system to unequivocally identify samples between sequence files and tables, and different naming systems are mixed. The first step is to fix that problem and use the same sample naming system in the full project. 

## Envision dataset.

### Samples identification and names.

In the case of the Envision project, the fastq files in the ENA servers under the same project are a mix of different samples and includes 16S and 18S rRNA genes from the field but also from incubation samples. As we are only interested in the 16S fastq files of field samples, we need to separate those from the other files prior to any analysis. 

To identify and separate the field samples of 16S rRNA gene from other kind of samples, we can use the information of the sequencing form that was an Excel file. After manually selectiong the samples of interest, create a list: `sample_list_raw_env.tsv`. 

However, the sample names of the sequencing form are different to the sample names used in the fastq files, sequence files information table and environmental metadata table. To fix this issue, first, we have to create a table of correspondences between sample names of interest with sample identification on environmental metadata and fastq files. The environmental metadata table contains sample IDs with the format: `ENV_Prok_number` that it is used as template to create sampleIDs for both Envision and Dimension data, but removing the "_" characters as they are problematic for downstream analysis.

```{r}
# Import list of sample names and labels for field samples of sterivex filters 
# amplified with 16S primers

# List produced from the sequencing form
raw_list <- read.table("data/files_env/sample_list_raw_env.tsv", sep="\t",header=TRUE)

# Change sample name to match those in fastq file tables ENV_Prok_###
library(tidyr)
lista <- raw_list %>%
  separate(sample_name, sep = "_", into = c("gene","number"))
# change number format to three digits
lista$number <- sprintf("%03d",as.integer(lista$number))

lista$sample_alias <- paste("ENV_Prok_", lista$number, sep="")

# Extract information of station, month, day and depth from labels
# to later merge with metadata
lista$sample_label_2 <- gsub("__","_", lista$sample_label) 
# labels from February do not follow the same order than the rest of the campaigns. 
# change them separately
lista_feb <- lista[1:44,]
lista_apraug <- lista[45:131,]

lista_feb <- lista_feb %>%
  separate(sample_label_2, sep = "_", into = c("gen2","month","st","day","prof_label"))
lista_apraug <- lista_apraug %>%
  separate(sample_label_2, sep = "_", into = c("gen2","month","day","st","prof_label"))

lista_feb <- lista_feb[,names(lista_apraug)]
lista <- rbind(lista_feb,lista_apraug)

# change variables to match metadata coding
# library(car)
source("./analysis/scripts/recode.r")

lista$month <- recode(lista$month, "'1602' = 'Feb';
                                       '1604' = 'Apr';
                                       '1608' = 'Aug'")

lista$day <- recode(lista$day, "'D1' = '1';
                                     'D3' = '3';
                                     'D5' = '5';
                                     'D7' = '7'", as.numeric = TRUE)

lista$st <- recode(lista$st, "'st3' = '3';
                                   'st6' = '6'", as.numeric = TRUE) 

# create a new variable with sample names whitout underscores (Mothur format)
lista$sample_name <- gsub("_","",lista$sample_alias)

# clean up data table
lista <- lista[,c("sample_alias","sample_name","month","day","st","prof_label")]
# Add sample that was not sequenced: August - st6 - D1, P4 
lista[nrow(lista)+1,] <- c("not_sequenced","notsequenced","Aug",7,6,"P4")

write.table(lista, "analysis/intermediate/sample_list_env.tsv", 
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")

rm("lista_apraug", "lista_feb", "raw_list")

###############################################
# Create table of correspondences between 
# sampleIDs and old fastq filenames
###############################################

# load table with project information from public databases - contains more samples than needed
sra <- read.table("data/files_env/filereport_read_run_PRJEB36188.tsv", header =TRUE, sep="\t")

# There is one sample less in the SRA filtered file because the list have all days with metadata
# but there is one not sequenced (August, st6, D7 P4). This sample is already included in "lista"
mysra <- subset(sra, sra$sample_alias %in% lista$sample_alias)

# drop unneeded columns
cols <- c("sample_alias","sample_accession","experiment_accession","run_accession","fastq_ftp","submitted_ftp")
mysra <- mysra[,cols]

# separate fastq read1 and read2 file names
mysra <- mysra %>%
  separate(fastq_ftp, sep = ";", into = c("ena_fastq_R1","ena_fastq_R2"))
mysra <- mysra %>%
  separate(submitted_ftp, sep = ";", into = c("original_fastq_R1","original_fastq_R2"))

# remove paths
mysra$original_fastq_R1 <- gsub("ftp.sra.ebi.ac.uk/vol1/run/ERR\\d+/ERR\\d+/","", mysra$original_fastq_R1) 
mysra$original_fastq_R2 <- gsub("ftp.sra.ebi.ac.uk/vol1/run/ERR\\d+/ERR\\d+/","", mysra$original_fastq_R2) 
mysra$ena_fastq_R1 <- gsub("ftp.sra.ebi.ac.uk/vol1/fastq/ERR\\d+/\\d+/ERR\\d+/","",mysra$ena_fastq_R1)
mysra$ena_fastq_R2 <- gsub("ftp.sra.ebi.ac.uk/vol1/fastq/ERR\\d+/\\d+/ERR\\d+/","",mysra$ena_fastq_R2)

# save the overall table with fastq files per sample, ENA and original
write.table(mysra, "analysis/intermediate/list_samples_fastq_files_original_and_ena_env.tsv", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# melt table to get one row per ENA fastq file to later use for selection and renaming
mysra_ena <- reshape2::melt(mysra[,-c(7:8)], value.name = "file", variable.name = "read",
                            id.vars = c("sample_alias","sample_accession","experiment_accession","run_accession"))

# create new fastq names with mothur format and sample ids
mysra_ena$read <- gsub("ena_fastq","", mysra_ena$read)
mysra_ena$sample <- gsub("_","",mysra_ena$sample_alias)
mysra_ena$newname <- paste(mysra_ena$sample, mysra_ena$read, sep ="")
mysra_ena <- mysra_ena[,c("sample", "sample_alias","sample_accession","experiment_accession","run_accession","file","newname")]

write.table(mysra_ena, "analysis/intermediate/list_samples_fastq_files_ena_only_env.tsv", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

```

This new file has this structure:
```{bash}
sample	sample_alias	sample_accession	experiment_accession	run_accession	file	newname
ENVProk001	ENV_Prok_001	SAMEA6473161	ERX3814636	ERR3813281	ERR3813281_1.fastq.gz	ENVProk001_R1
ENVProk002	ENV_Prok_002	SAMEA6473162	ERX3814637	ERR3813282	ERR3813282_1.fastq.gz	ENVProk002_R1
ENVProk003	ENV_Prok_003	SAMEA6473163	ERX3814638	ERR3813283	ERR3813283_1.fastq.gz	ENVProk003_R1
ENVProk004	ENV_Prok_004	SAMEA6473164	ERX3814639	ERR3813284	ERR3813284_1.fastq.gz	ENVProk004_R1
ENVProk005	ENV_Prok_005	SAMEA6473165	ERX3814640	ERR3813285	ERR3813285_1.fastq.gz	ENVProk005_R1
ENVProk006	ENV_Prok_006	SAMEA6473166	ERX3814641	ERR3813286	ERR3813286_1.fastq.gz	ENVProk006_R1
```


### Environmental metadata pre-process

Now that we had this table with name correspondences, the next step was to correct the environmental metadata table. 

In the case of Envision project, the environmental metadata consisted on two kind of tables. On one hand, three csv files further called raw metadata (`Metadata_ENV1602.csv`, `Metadata_ENV1604.csv` and, `Metadata_ENV1608.csv`) that included a set of measurements taken in February, April and August respectively (i.e. nitrate and nitrite concentrations or different cell counts). These files do not have the same variable or sample names, not even the same structure. On the other hand, there is a set of tables with data [from a CTD][https://en.wikipedia.org/wiki/CTD_(instrument)], that is an oceanographic instrument that measures conductivity (C), Temperature (T) and Depth (D) among other variables. In this case, I have six tables, corresponding to the three different months and two sampling places. This data needs to be processed and joined with the raw metadata to produce a single table with all environmental variables for all samples with their corrected names.

The script following corrects the table's structure, naming systems and join the different tables into a single one. 

Prior to load in R, CTD files were saved as tab-separated text files and encoded to UTF8.

```{r}
###########################################################
#
# Envision metadata correction 
# 
###########################################################

setwd("/data/mcm/nfernandez/envisdim")

# Load data and preprocess 
###########################################################
# Load CDT data adding st and month when needed ####
fctd <- c("analysis/intermediate/CTD_data/CTD_Feb_ST3.utf8.txt", 
          "analysis/intermediate/CTD_data/CTD_Feb_ST6.utf8.txt",
          "analysis/intermediate/CTD_data/CTD_Apr_ST3.utf8.txt",
          "analysis/intermediate/CTD_data/CTD_Apr_ST6.utf8.txt",
          "analysis/intermediate/CTD_data/CTD_Aug_ST3.utf8.txt",
          "analysis/intermediate/CTD_data/CTD_Aug_ST6.utf8.txt")

ctd <- lapply(fctd, read.table, sep = "\t", header = TRUE, dec = ".")
# In February tables variable names presion and depth are twisted, 
# the order is depth - presion, as in April and August tables,
# but it is no a problem as this will be corrected later.

# remove extra variable, Long, in April and August data
ctd[3:6] <- lapply(ctd[3:6],`[`,-3)

# Add station variable to February data
ctd[[1]]$st <- rep(3, nrow(ctd[[1]]))
ctd[[2]]$st <- rep(6, nrow(ctd[[2]]))
# put first to match the other datasets
ctd[1:2] <- lapply(ctd[1:2],`[`, c(11,1:10))

# Add month variable
for (i in 1:length(ctd)) {
  if (i %in% 1:2) {
    ctd[[i]]$month <- rep(c("Feb"), nrow(ctd[[i]]))
  } else if (i %in% 3:4) {
    ctd[[i]]$month <- rep(c("Apr"), nrow(ctd[[i]]))
  } else {
    ctd[[i]]$month <- rep(c("Aug"), nrow(ctd[[i]]))
  }
}
# put it first
ctd <- lapply(ctd,`[`, c(12,1:11))

# Variable names have different formats, homogenize them 
# (this step takes care of presion - depth problem in February).
ctd <- lapply(ctd, setNames, names(ctd[[3]]))
# Join all ctd tables into a single one
ctdall <- do.call("rbind", ctd)

# remove data from days not sampled
ctdall <- subset(ctdall, subset = Day %in% c(1,3,5,7))

# Change detph from negative to positive numbers
ctdall$depth <- ctdall$depth * -1

rm(ctd, fctd,i)

# Load raw metadata from excel files saved as csv ####
rawFeb <- read.table("data/files_env/Metadata_ENV1602.csv", header = TRUE, sep = "\t")
rawApr <- read.table("data/files_env/Metadata_ENV1604.csv", header = TRUE, sep = "\t")
rawAug <- read.table("data/files_env/Metadata_ENV1608.csv", header = TRUE, sep = "\t")

# reorder August variables, they do not match the other campaigns
rawAug <- rawAug[,c(1:8,11,9,10,12:35)]

# change variables names to more R friendly versions
nombres <- c("day","st","muestras","collection","prof","no3","no2","nh4","TIN",
             "po4","sio2","DOC","TN","DON","chlaLow3","chlaUpp3","chlaTotal", 
             "AB","BB","AFSyne","AFProc","AFsp","AFlp","BFSyne", "BFProc",
             "BBsp","BBlp","PPupp3","eePPup3","PPlow3","eePPup3","PPtotal",
             "eePPtotal","BP","eeBP") 

names(rawFeb) <- nombres
names(rawApr) <- nombres
names(rawAug) <- nombres

# add month variable
rawFeb$month <- rep(c("Feb"), nrow(rawFeb))
rawApr$month <- rep(c("Apr"), nrow(rawApr))
rawAug$month <- rep(c("Aug"), nrow(rawAug))

# Merge the three campaigns and reorder variables
rawdata <- rbind(rawFeb, rawApr, rawAug)
rawdata <- rawdata[,c(36,1:35)]

# Remove stations other than 3 and 6 and days 1,3,5,7
rawdata <- subset(rawdata, subset = st %in% c(3,6) & day %in% c(1,3,5,7))
row.names(rawdata) <- c(1:nrow(rawdata))

rm("rawApr","rawAug","rawFeb", "nombres")

# Add CTD data to rawdata
###########################################################

ctdsel <- list()

for (i in 1:nrow(rawdata)){
  # select ctd data for that campaign, station and day only
  ctd_part <- subset(ctdall, month == rawdata$month[i] & st == rawdata$st[i] & Day == rawdata$day[i]) 
  # get the depth of the DNA sample
  profun <- rawdata$prof[i]
  # find the index of the closest depth in ctd data subset
  j <- which(abs(ctd_part$depth - profun) == min(abs(ctd_part$depth - profun)))
  # select the row for the index
  ctdsel[[i]] <- ctd_part[j,-c(1,2,3,5)]
  
  rm("i","j","ctd_part","profun")
}

# create a dataframe for selected ctd data
ctdsel <- do.call("rbind", ctdsel)
row.names(ctdsel) <- 1:nrow(ctdsel)
# keep ctd real depth, it is not the same then dna data in all cases
names(ctdsel)[1]<-"ctdProf"

#merge with the rest of metadata 
metadata <- cbind(rawdata, ctdsel)

# Add prof_label variable to merge with 16S sample list
# separate muestras by profundity labels
etiquetas <- c()
etiquetas[grep("prof 1", metadata$muestras)] <- "P1"
etiquetas[grep("p1", metadata$muestras)] <- "P1"
etiquetas[grep("prof 2", metadata$muestras)] <- "P2"
etiquetas[grep("p2", metadata$muestras)] <- "P2"
etiquetas[grep("prof 3", metadata$muestras)] <- "P3"
etiquetas[grep("p3", metadata$muestras)] <- "P3"
etiquetas[grep("prof 4", metadata$muestras)] <- "P4"
etiquetas[grep("p4", metadata$muestras)] <- "P4"
etiquetas[grep("prof 5", metadata$muestras)] <- "P5"
etiquetas[grep("p5", metadata$muestras)] <- "P5"
etiquetas[grep("prof 6", metadata$muestras)] <- "P6"
etiquetas[grep("p6", metadata$muestras)] <- "P6"
etiquetas[grep("prof 7", metadata$muestras)] <- "P7"
etiquetas[grep("p7", metadata$muestras)] <- "P7"
metadata$prof_label <- etiquetas

# Add names from sequenced samples
lista <- read.table("analysis/intermediate/sample_list_env.tsv", header = TRUE, sep = "\t")
metadata <- merge(metadata, lista, by = c("month","day","st","prof_label"))

# reorder a bit
metadata <- metadata[, c(46,47,1:6,38,7:37,39:45)]

# save to file
write.table(metadata, file = "analysis/intermediate/metadata_corrected_full_table_env.tsv",
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")
write.csv2(metadata, file = "analysis/intermediate/metadata_corrected_full_table_env.csv")

# Two samples from August st6 does not have 16S Data, remove them and save data again
# samples= d1, prof1 = 5m and d7 prof4  -> names shared between April and August
filas <- as.numeric(rownames(subset(metadata, month == "Aug" & muestras %in% c("st6 day1 p1_ 19", "st6 day7 p4_90"))))
metadata_16S <- metadata[-filas,]

write.table(metadata_16S, file = "analysis/intermediate/metadata_corrected_16S_table_env.tsv",
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")
write.csv2(metadata_16S, file = "analysis/intermediate/metadata_corrected_16S_table_env.csv")
```


### Rename fastq files

Using the correspondence names file previously created for Envision bioproject, select and rename the fastq files.

```{bash}
#!/bin/bash

# rename (new sample ID) and copy 16S fastq files  
# Envision project

# Define pathways and file names corresponding to 16S field samples
dir1='data/files_env/raw_fastq_env'
dir2='analysis/intermediate/raw_16_fastq'
datain='analysis/intermediate/list_samples_fastq_files_ena_only_env.tsv'
table='analysis/intermediate/list_samples_fastq_files_ena_only_env_no_head.tsv'

# Remove header
cat $datain | sed '1d' > $table

# Rename files while selecting them
while read -r sample sample_alias sample_accession experiment_accession run_accession file newname ; do
   cp $dir1/$file $dir2/$newname.fastq.gz
done < $table
```



## Dimension dataset

### Sample identification, names and environmental metadata preprocess

This project contains data from two different years. The naming systems used for the fastq files on each year were completely different. All environmental data was included in a single file, but this table had no sample IDs. The only way to identify the correspondence between fastq files and environmental data was through the sampling date, that was included in the fastq file name for the samples of the second year. For the samples of the first year, we had a correspondence table between a numeric identification included in the fastq name and the sampling date.

The script following script corrects these differences and creates a sample naming system similar to the one used in Envision project (**DIM_Prok_#**), while keeping the numeric identification of the first year samples. It also modifies the environmental metadata to include the new sample IDs.

```{r}
###########################################################
#
# Dimension - file names and metadata corrections
#
###########################################################

setwd("/data/mcm/nfernandez/envisdim")
# setwd("D:/OneDrive/TRABAJO.CSIC/prj.tfm")

# load fastq filenames
raw_list <- read.table("./data/files_dim/filenames_16S_dimension.txt", 
                        header = FALSE)
names(raw_list) <- "filename"

# Sample correspondence list for 2014
corresp <- read.csv2("./data/files_dim/DIMENSION_sample_correspondence.csv")
# library(car)
source("./analysis/scripts/recode.r")

corresp$month <- recode(corresp$month, "'enero' = 'Jan';
                                             'feb' = 'Feb';
                                             'marzo' = 'Mar';
                                             'abril' = 'Apr';
                                             'mayo' = 'May';
                                             'junio' = 'Jun';
                                             'sept' = 'Sep';
                                             'oct' = 'Oct'")
corresp$depth <- recode(corresp$depth, "'0' = '5'", as.numeric = TRUE)

# Load metadata and reorganize
rawmetad <- read.csv2("./data/files_dim/metadatosDIMENSION.csv")
rawmetad$month <- gsub("\\d+.","",rawmetad$month)
rawmetad$year <- as.numeric(gsub("\\d+/","", rawmetad$Dat))
rawmetad$Dep <- as.numeric(rawmetad$Dep)
rawmetad <- rawmetad[,c(1,52,2:51)]

# Get info from filenames
raw_list$filename <- gsub("raw_fastq_dim/","",raw_list$filename)

# filenames from 2014 and 2015 have different formats, divide the list
list14 <- data.frame(filename = raw_list[1:34,])
list15 <- data.frame(filename = raw_list[35:74,])
# add year
list14$year <- rep(2014, nrow(list14))
list15$year <- rep(2015, nrow(list15))

library(tidyr)
# Divide filenames and drop unneeded information
list14 <- list14 %>%
  separate(filename, sep = "_", into = c("info","reads"), remove = FALSE) %>%
  separate(info, sep="-", into = c("num1","sample","prim1","prim2"))
list14 <- list14[,c("year","sample","reads", "filename")]

# Add sample information from correspondence table
list14$month <- c(NA)
list14$depth <- c(NA)

for (i in 1:nrow(list14)){
  for (j in 1:nrow(corresp)){
    if (list14$sample[i] == corresp$sample[j]){
      list14$month[i] <- corresp$month[j]
      list14$depth[i] <- corresp$depth[j]
    }
  }
  
}
# Reorganize 
list14 <- list14[,c("year", "month","depth","reads", "filename", "sample")]

# Repeat process for year 15, adapting the steps to data different structure
list15 <- list15 %>%
  separate(filename, sep = "_", into = c("info","reads"), remove = FALSE) %>% 
  separate(info,"-",into = c("num1","bac","month","depth","pico","prim1","prim2"))
list15 <- list15[,c("year", "month","depth","reads", "filename")]

list15$month <- recode(list15$month,"'abr15' = 'Apr'; 
                                          'ago15' = 'Aug';
                                          'feb15' = 'Feb';
                                          'jul15' = 'Jul';
                                          'jun15' = 'Jun';
                                          'mar15' = 'Mar';
                                          'may15' = 'May';
                                          'nov15' = 'Nov';
                                          'oct15' = 'Oct';
                                          'sep15' = 'Sep'")
list15$depth <- recode(list15$depth,"'0m' = '5'; '30m' = '30'", as.numeric=TRUE)

list15$sample <- c(NA)

for(i in 1:nrow(list15)){
  for(j in 1:nrow(rawmetad)){
    if(list15$year[i] == rawmetad$year[j] && list15$month[i] == rawmetad$month[j] && list15$depth[i] == rawmetad$Dep[j]){
      list15$sample[i] <- rawmetad$sample[j]
    }
  }
}

# Merge lists of 2014 and 2015
lista <- rbind(list14, list15)
lista$reads <- gsub(".fastq","", lista$reads)
lista <- lista[,c("sample", "year","month","depth","reads","filename")]

# Create new sample IDs 
lista$newname <- paste("DIMprok",as.character(lista$sample), "_", lista$reads, sep = "")

# Save file to rename fastq files
write.table(lista, 
            "./analysis/intermediate/list_samples_fastq_files_dim.tsv", 
            sep="\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

# Add new sampleIDs to environmental metadata and save
metadata <- rawmetad
metadata$sample <- paste("DIMProk", as.character(rawmetad$sample), sep = "")

# recode Season variable
metadata$Sea <- recode(metadata$Sea, "'Otono' = 'Fall';
                                      'Invierno' = 'Winter';
                                      'Primavera' = 'Spring';
                                      c('verano','Verano') = 'Summer'")
# Save data into a file
write.table(metadata, 
            "./analysis/intermediate/metadata_corrected_dim.tsv", 
            sep="\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)
```


This script generated two tables. 
The file `list_samples_fastq_files_dim.tsv` that will be use to rename the fastq files with the new sample IDs:
```{bash}
sample	year	month	depth	reads	filename	newname
10	2014	2	30	R1	4316-10-515yF-926R_R1.fastq	DIMprok10_R1
10	2014	2	30	R2	4316-10-515yF-926R_R2.fastq	DIMprok10_R2
13	2014	5	5	  R1	4316-13-515yF-926R_R1.fastq	DIMprok13_R1
13	2014	5	5	  R2	4316-13-515yF-926R_R2.fastq	DIMprok13_R2
1	  2014	3	5	  R1	4316-1-515yF-926R_R1.fastq	DIMprok1_R1
1	  2014	3	5	  R2	4316-1-515yF-926R_R2.fastq	DIMprok1_R2
```

And the corrected environmental metadata file `metadata_corrected_dim.tsv`.


### Rename fastq files

As done with Envision metadata, a bash script is used to rename the fastq files.

```{bash}
#!/bin/bash

# rename (new sample ID) and copy 16S fastq files  
# Dimension project

# define paths
dir1='data/files_dim/raw_fastq_dim'
dir2='analysis/intermediate/raw_16_fastq'
datain='analysis/intermediate/list_samples_fastq_files_dim.tsv'
table='analysis/intermediate/list_samples_fastq_files_dim_no_head.tsv'

# remove headers
cat $datain | sed '1d' > $table

# Rename files while selecting them
while read -r sample year month depth reads filename newname ; do
   cp $dir1/$filename $dir2/$newname.fastq
done < $table

# remove duplicated sample
rm $dir2/DIM*49*fastq
```

             
             
  
             
             
             
             