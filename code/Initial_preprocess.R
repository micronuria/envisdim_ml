
# Author: Nuria Fernández-Gonzalez
# Date: 29-05-2023

###########################################
#
# INITIAL DATA PREPROCESS
#
###########################################

#---------------------------------------------
# Description:
# This script:
#
#  1 - reads the raw ASVs calculated with DADA2, corrects
#      sample names from Dimension campaign and removes a problematic
#      sample from Envision campaign for which the sequencing did not work. 
#
#  2 - Loads and joins metadata from Envision and Dimension campaigns.
#
#  3 - Loads chlorophyll based sample groups.
# 
#  4 - Adds campaign, event-campaign and season categorical variables.
#---------------------------------------------

###########################################
# Load dependencies
###########################################

library(tidyverse)


###########################################
# Load ASVs count raw and taxonomy data
###########################################

# 1 - Correct sample names of Dimension project to meet metadata and Envision
# sample formats.
# 2 - Delete ENVISION sample with sequencing problems:
# ENVProk141 - Ocean, August, Day 7, depth 1
# 3 - Table to long format
# 4 - Load taxonomy

all_asv <- read_tsv("raw_data/asvs_raw_counts.tsv",
            col_types = cols(sample = col_character(),
                             .default = col_integer())) %>%
  mutate(sample = str_replace_all(sample, "DIMprok", "DIMProk")) %>%
  filter(!(sample == "ENVProk141"))%>% 
  pivot_longer(-sample, names_to = "asv", values_to = "count")

# load taxonomy 
tax_raw <- read_tsv("raw_data/asvs_taxonomy.tsv")

###########################################
# Load metadata for Envision and Dimension
###########################################

# remove variables not needed for the initial exploratory analysis
met_raw_env <- read_tsv("raw_data/metadata_corrected_16S_table_env.tsv") %>%
  mutate(year = 2016) %>%
  rename(sample = sample_name, depth = prof) %>%
  select(sample, st, year, month, depth)

met_raw_dim <- read_tsv("raw_data/metadata_corrected_dim.tsv") %>%
  rename(depth = Dep) %>%
  mutate(st = 3) %>%  # Add station data
  select(sample, st, year, month, depth)

# join both metadata tables and add season and 
# depth level categorical variables
met_raw <- bind_rows(met_raw_dim, met_raw_env) %>%
  mutate(season = case_when(month %in% c('Jan', 'Feb','Mar') ~'Winter',
                            month %in% c('Apr','May','Jun') ~ 'Spring',
                            month %in% c('Jul','Aug','Sep') ~'Summer',
                            TRUE ~ 'Fall')) %>%
  mutate(depth_lev = case_when(depth %in% c('5','8','10','15','20','25','30') ~ 'surface',
                               depth %in% c('40','50','55') ~ 'intermediate',
                               depth %in% c('60', '70', '75', '80', '85', 
                                            '90', '100', '110', '125', '130',
                                            '150', '160', '200') ~ 'deep'))

###########################################
# Load Chlorophyll groups
###########################################

# Read sample grouping according to bloom events
# and drop chlorophyll continuous data
chlclass <- read_tsv("raw_data/chlclass.tsv") %>%
  select(-chl)

# check number of cases on bloom and no-bloom categories
# Neet to keep this in mind to train the ML algorithm
# as the data is imbalanced
chlclass %>% count(event)
# # A tibble: 2 × 2
# eve nt      n
# <chr>  <int>
# 1 bloom     31
# 2 normal   135

# Creates a campaign variable and a 
# campaign x event (outcomes) variable
chlclass <- chlclass %>%
  mutate(
    campaign = case_when(
      str_detect(sample, 'ENV')  ~ 'envision',
      TRUE ~ 'dimension')) %>%
    mutate(
    gr_event = case_when(
      str_detect(sample, 'ENV') & event == 'bloom' ~ 'bloom_env',
      str_detect(sample, 'ENV') & event == 'normal' ~ 'normal_env',
      str_detect(sample, 'DIM') & event == 'bloom' ~ 'bloom_dim',
      TRUE ~ 'normal_dim'))

# Check number of cases on normal and bloom categories 
# by campaign 
chlclass %>%
  with(table(gr_event))
# gr_event
# bloom_dim  bloom_env normal_dim normal_env
# 10         21         26        109

# Check number of cases by campaign 
chlclass %>%
  with(table(campaign))
# campaign
# dimension  envision 
# 36       130 


###########################################

# Clean environment
rm(list = c("met_raw_dim", "met_raw_env"))
# Save objects
save(list = ls(), file = "processed_data/initial_preprocess.RData")


