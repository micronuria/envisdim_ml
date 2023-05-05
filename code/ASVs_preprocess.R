###########################
# DATA PREPROCESS
##########################


library(tidyverse)


# load ASVs count raw data
read_tsv("raw_data/asvs_raw_counts.tsv", 
         col_types = cols(sample = col_character(),
                          .default = col_integer())) 

%>% 
  pivot_longer()



### create a dataframe with taxonomy


# add taxonomy columns
asvs_raw_t <- cbind(taxa.sp.print, asvs_raw)
#asvs_raw_t$tax <- with(asvs_raw_t, paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "__"))


# check NAs in taxonomy higest level
table(is.na(asvs_raw_t$Kingdom))
# FALSE  TRUE 
# 8178     6 

table(asvs_raw_t$Kingdom)
# Archaea  Bacteria Eukaryota 
#      665      7508         5 

dim(asvs_raw_t)
# [1] 8184  174

# Remove Kingdom: "Eukaryota" and NA
asvs_raw_t <- subset(asvs_raw_t, !(Kingdom == "Eukaryota"))
dim(asvs_raw_t)
# [1] 8173  174
table(is.na(asvs_raw_t$Kingdom))
# FALSE 
# 8173 

### TEST - grouping ASVs by tax
# to long format
asvs_raw_t.l <- melt(asvs_raw_t, id.vars = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                     value.name = "counts",
                     variable.name = "sample")

# group by taxonomy and sample, sum counts
asvs_grouped.l <- asvs_raw_t.l %>% 
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, sample) %>%
  summarise(value = sum(counts))
# to wide format
asvs_test.w <- asvs_grouped.l %>%
  pivot_wider(names_from = sample, values_from = value)

asvs_test.w$tax <- with(asvs_test.w, paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "__"))

asvs_test <- t(asvs_test.w[,c(175,8:174)])

# Add variable with sample classification
chlclass <- read.table("/data/mcm/nfernandez/envisdim/analysis/intermediate/env_metadata/chlclass.tsv",
                       sep = "\t", header = TRUE)


