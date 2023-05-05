############################################
# Write data tables from DADA2 results
############################################

# load DADA2 results
load("raw_data/DADA2_CR30.No_primers.RData")

# transpose count data matrix
asvs_raw <- as.data.frame(seqtab.nochim)

# separate sequences from count data
asvs_seqs <- data.frame(seq = names(asvs_raw))
# create ASVs names
names(asvs_raw) <- paste("asv",sprintf("%04d", 1:ncol(asvs_raw)), sep="")
asvs_seqs$asv <-names(asvs_raw)
asvs_seqs <- asvs_seqs[,c("asv","seq")]

# sample names to a column
asvs_raw$sample <- rownames(asvs_raw)
asvs_raw <- asvs_raw[, c(ncol(asvs_raw), 2:ncol(asvs_raw)-1)]

# taxonomy
taxonomy <- as.data.frame(taxa.sp.print)
taxonomy$asv <- names(asvs_raw)[2:ncol(asvs_raw)]

# save to files
# counts
write.table(asvs_raw, file = "raw_data/asvs_raw_counts.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) 
# sequences
write.table(asvs_seqs, file = "raw_data/asvs_seqs.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# taxonomy
write.table(taxonomy, file = "raw_data/asvs_taxonomy.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
