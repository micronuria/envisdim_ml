############################################
# Write data tables from DADA2 results
############################################

# load DADA2 results
load("raw_data/DADA2_CR30.No_primers.RData")

# transpose count data matrix
asvs_raw <- as.data.frame(t(seqtab.nochim))

# separate sequences from count data
asvs_seqs <- rownames(asvs_raw)
rownames(asvs_raw) <- 1:nrow(asvs_raw)

# save to files
# counts
write.table(asvs_raw, file = "raw_data/asvs_raw_counts.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) 
# sequences
write.table(asvs_seqs, file = "raw_data/asvs_seqs.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# taxonomy
write.table(taxa.sp.print, file = "raw_data/asvs_taxonomy.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
