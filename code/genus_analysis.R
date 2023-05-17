source('code/ASVs_preprocess.R')

library('mikropml')

# drop variables not needed to train the model
chl_genus_data <- composite_wide %>% 
  select(-sample, -chl)


