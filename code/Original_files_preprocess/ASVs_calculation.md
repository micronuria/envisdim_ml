
# ASVs calculation

Calculation of Amplicon Sequence Variants (ASVs) with DADA2 in R.


## DADA2 installation

First, we create a conda environment to install and R version compatible with latest BiocManager and install both.
```{bash}
conda create -n dada2
conda activate dada2
conda install -c conda-forge r-base=4.2.0 # compatible with latests BiocManager stable version
```

After launching R, install dada2 package
```{r}
install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("dada2", version = "3.16")
```


## ASVs calculation

Next, we can calculate the ASVs using the preprocessed sequences.
Script based on my coworker Diego Jimenez's code and on [DADA2 tutorial](https://benjjneb.github.io/dada2/index.html)).

```{r}
library(dada2)

setwd("/media/disk5/nfernandez/envisdim_trim/fastq_for_asvs")


path.cr30 <- '/media/disk5/nfernandez/envisdim_trim/fastq_for_asvs'
list.files(path.cr30)
print(path.cr30)

# create vector of files
cr30 <- sort(list.files(path.cr30, pattern=".fastq", full.names = FALSE))

# Obtain sample name and associate them to file names
sample.names <- sapply(strsplit(basename(cr30), "primerfree_"), `[`, 2)
sample.names <- sapply(strsplit(sample.names, ".preR.fastq"), `[`, 1)
names(cr30) <- sample.names
names(cr30)

error.cr30 <- learnErrors(cr30, qualityType="FastqQuality", 
                          multithread=TRUE, randomize=TRUE)

png(file = "plot_quality.no_primers_DADA2_CR30.FastqQuality.png", 
    bg = "transparent", width=1200,height=684)
plotErrors(error.cr30, nominalQ=TRUE)
dev.off()

print("Learn Error rates have been finished")

# Checkpoint 1
save.image(file="DADA2_CR30.No_primers.RData")
print("Checkpoint 1 achieved")

dada.cr30 <- dada(derepFastq(cr30, qualityType="FastqQuality"), 
                  err=error.cr30, pool=TRUE, multithread=TRUE)
print("The dada object has been created")

# Checkpoint 2
save.image(file="DADA2_CR30.No_primers.RData")
print("Checkpoint 2 achieved")

# Construction of the table for further steps
seqtab <- makeSequenceTable(dada.cr30)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

print("The percentage of sequences that has been passed the chimera removing is")
sum(seqtab.nochim)/sum(seqtab)

# NÂº ASVs Recovered/sample
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dada.cr30, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "nonchim")
rownames(track) <- sample.names
track

# Checkpoint 3
save.image(file="DADA2_CR30.No_primers.RData")
print("Checkpoint 3 achieved")

# Taxonomic assignment - SILVA v.138 (Domain - Phylum - Order - Class - Family - Genus)
taxa <- assignTaxonomy(seqtab.nochim, "/media/disk5/nfernandez/envisdim_trim/SILVAv132_DADA2/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
print("Taxonomy assigned (until Genus level)")

# Checkpoint 4
save.image(file="DADA2_CR30.No_primers.RData")
print("Checkpoint 4 achieved. You have a taxonomy assignment for your ASVs in taxa object")

# Getting specific assignment for taxonomy
taxa.sp <- addSpecies(taxa, "/media/disk5/nfernandez/envisdim_trim/SILVAv132_DADA2/silva_species_assignment_v138.1.fa.gz")

# Checkpoint 5
save.image(file="DADA2_CR30.No_primers.RData")
print("Checkpoint 5 achieved. The Specie level has been added to your ASV table in taxa.sp object")

# Removing sequence rownames for display only
taxa.sp.print <- taxa.sp # Removing sequence rownames for display only
rownames(taxa.sp.print) <- NULL
taxa.sp.print

# Final Checkpoint 
save.image(file="DADA2_CR30.No_primers.RData")
print("Final Checkpoint achieved. The DADA2 Rscript has finished.")
print("From this point, you could export your tables and work with them with VEGAN, PHYLOSEQ or any other you want.")


getwd()
```
