###########################
# DATA PREPROCESS
##########################


library(tidyverse)

set.seed(20210401)

# load ASVs count raw data
# There is an ENVISION sample with sequencing problems
# ENVProk141 - Ocean, August, Day 7, depth 1
# remove it

all_asv <- read_tsv("raw_data/asvs_raw_counts.tsv", 
         col_types = cols(sample = col_character(),
                          .default = col_integer())) %>% 
         mutate(sample = str_replace_all(sample, "DIMprok", "DIMProk")) %>% # correct an error in sample names of Dimension project
  filter(!(sample=="ENVProk141")) %>% # remove problematic sample
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

# read sample grouping according to bloom events
chlclass <- read_tsv("raw_data/chlclass.tsv") 

# check number of cases on bloom and no-bloom categories
# Neet to keep this in mind to train the ML algorithm
# as the data is imbalanced

chlclass %>% count(event)
# # A tibble: 2 × 2
# eve nt      n
# <chr>  <int>
# 1 bloom     31
# 2 normal   135

# join count and taxonomy and ...
composite <- inner_join(all_asv, taxonomy, by = "asv") %>%
  group_by(sample, taxonomy) %>%                        # group by genus
  summarize(count = sum(count), .groups="drop") %>% 
  group_by(sample) %>%                                  # calculate relative abundances
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         # get rid of grouping structure on data
  select(-count) %>%                                    # get rid of count column
  inner_join(., chlclass, by="sample")                  # add data about chlorophyll groups



# to wide format
composite_wide <- composite %>%
  pivot_wider(names_from = taxonomy, values_from = rel_abun)






