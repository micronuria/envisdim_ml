# Repository of TFM
## TITLE

Code developed during the TFM can be found in the code folder that contains:

#### Original_files_preprocess:
- Data dowload from SRA:
    - `Data_collection.md`
    - `sra_explorer_download.sh`
- Environmental metadata and sequence files preparation:
    - `Metadata_files_preprocess.md`
- Raw sequences preprocess:
    - `Sequences_preprocess.md`
- ASVs calculation:
    - `ASVs_calculation.md`

#### Chlorophyll_groups
- PCAs of environmental variables: 
    - `PCA_environmental_data.R`
- Data distribution according to different percentiles of total chlorophyll concentrations:
    - `Percentiles.R`

#### Main folder scripts:
- Data preparation:
    - Preparation of DADA2 results: `DADA2_tables.R`
    - Merge of Envision and Dimension environmental metadata and sample labelling: `Initial_preprocess.R`
- Data exploration:
    - ASVs without filtering exploration: `ASVs_no_filter_exploration.R`
    - ASVs filtering and exploration: `ASVs_filter_exploration.R`
    - Clusters generation and exploratory analyses: `Clusters_grouping_exploration.R`
- Train, validation with 5 and 10 k-fold repeated cross-validation 
    and performance test of Random Forest models:
    - Clusters: `Clusters_model_1.R`
    - Filtered ASVs: `ASVs_model_1.R`
    - Unfiltered ASVs: `ASVs_no_filter.R`
    - Validation and performace metrics: `RF_models_results.R`

        