# Code Repository
## Application of Machine Learning Methods to Predict Phytoplankton Blooms and Determine Microbial Biomarkers Using Marine Microbiomes

Code developed during the TFM entitled *Application of Machine Learning Methods to Predict Phytoplankton Blooms and Determine Microbial Biomarkers Using Marine Microbiomes* can be found in the code folder that contains:

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

#### Chlorophyll_groups:
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
    and performance test of Random Forest and SVM models:
    - Random Forest with Clusters: `Clusters_model_1.R`
    - Random Forest with Filtered ASVs: `ASVs_model_1.R`
    - Random Forest with Unfiltered ASVs: `ASVs_no_filter_model_1.R`
    - SVM Radial with Filtered ASVs and validation results: `ASVs_SVM_model.r`
    - Validation and performance metrics for RF and SVM models: `RF_SVM_models_results.R`
    
- Synthetic data:
   - Data generation and models training: `ASVs_synthetic_RF.R`
   - Performance results of RF model with synthetic data: `ASVs_synthetic_RF_results.R`

- Features importance: `feature_importance.R`

        
