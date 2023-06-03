




# calculate relative abundances
# join count data, taxonomy, metadata, chlorophyll groups 
# and 
composite <- inner_join(all_asv, taxonomy, by = "asv") %>%
  group_by(sample, taxonomy) %>%                        # group by genus
  summarize(count = sum(count), .groups="drop") %>% 
  #group_by(sample) %>%                                  # calculate relative abundances
  #mutate(rel_abun = count / sum(count)) %>% 
  #ungroup() %>%                                         # get rid of grouping structure on data
  #select(-count) %>%                                    # get rid of count column
  inner_join(., chlclass, by="sample") %>%              # add data about chlorophyll groups
  inner_join(., met_raw)






# Preprocess genus clusters by removing near-zero variance 
# and highly correlated features, but not negatively correlated features
# for interpretability
clusters_prep <- composite_wide %>% 
  select(-sample, -gr_event) %>%
  mutate(campaign = as.factor(campaign),
         st = as.factor(st),
         year = as.factor(year),
         month = as.factor(month)) %>%
  preprocess_data(outcome_colname = "event",
                  method = c("center","scale"),              # already preprocessed
                  collapse_corr_feats = TRUE, 
                  group_neg_corr = FALSE,    
                  remove_var='nzv',          
                  prefilter_threshold = -1) # do not filter by prevalence




