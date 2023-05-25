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

# read sample grouping according to bloom events
chlclass <- read_tsv("raw_data/chlclass.tsv") %>% 
  select(-chl)

# join 
#asv_raw <- 
asv_relabun <- all_asv %>% # add data about chlorophyll groups
  group_by(sample) %>%                                  # calculate relative abundances
  mutate(rel_abun = count / sum(count)) %>% 
  ungroup() %>%                                         # get rid of grouping structure on data 
  select(-count) %>%                                    # get rid of count column
  inner_join(., chlclass, by="sample") %>%
  pivot_wider(names_from = asv, values_from = rel_abun)

library(mikropml)
library(purr)
library(tictoc)

get_results <- function(seed) {
  run_ml(asv_relabun,
         method = "rf",
         outcome_colname = "event",
         cv_times = 1,
         seed = seed)
}


# check time for a single run with all ASVs
tic()
iterative_results <- purrr::map(c(1,2,3), get_results)
toc()


library(future)
library(furrr)

future::plan(multisession)
tic()
iterative_results <- furrr::future_map(c(1,2,3), get_results)
toc()
