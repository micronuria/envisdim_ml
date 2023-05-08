###########################
# DATA PREPROCESS
##########################


library(tidyverse)


# load ASVs count raw data
all_asv <- read_tsv("raw_data/asvs_raw_counts.tsv", 
         col_types = cols(sample = col_character(),
                          .default = col_integer())) %>% 
  pivot_longer(-sample, names_to = "asv", values_to = "count")


# GROUPED DATA

# load taxonomy data up to Genus level only
tax_raw <- read_tsv("raw_data/asvs_taxonomy.tsv") %>%
  select(asv, Kingdom, Phylum, Class, Order, Family, Genus) 


# Group taxonomy to filter data
taxonomy <- tax_raw %>%
  unite("taxonomy", Kingdom:Genus, sep = ";", remove = TRUE) %>%
  mutate(taxonomy = str_replace(taxonomy, ";NA","_unclassified"),
         taxonomy = str_replace_all(taxonomy, ";NA", ""),
         taxonomy = str_replace_all(taxonomy, ".*;", ""))

# read metadata
chlclass <- read_tsv("raw_data/chlclass.tsv") 

# check number of cases on bloom and no-bloom categories
chlclass %>% count(event)
 

# join count and taxonomy and ...
composite <- inner_join(all_asv, taxonomy, by = "asv") %>%
  group_by(sample, taxonomy) %>%                        # group by genus
  summarize(count = sum(count), .groups="drop") %>% 
  group_by(sample) %>%                                  # calculate relative abundances
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         # get rid of grouping structure on data
  select(-count) %>%                                    # get rid of count column
  inner_join(., chlclass, by="sample")                  # add data about chlorophill groups

##### add metadata!!!


#####################################
# REMOVE eUK ANd ALL NA - chloroplasts
###################################

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

 \\(\  MP)(\d+)


print $1 $2