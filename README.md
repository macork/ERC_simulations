# Simulation_studies
Code for simulations for Michael Cork paper #1: Curving emissions: Comparing methods for evaluating exposure response curves

This respository contains code used to conduct simulation studies comparing methods for exposure-response curves. It also contains a data application to the Medicare data set. 

The main code scripts are as follows:
`0_launch.jobs.sh` is script to launch jobs on the FASSE cluster 

There are three main code directories in this repo:
1) Simulation studies which has the code 
2) Dan IPTW example which is a slightly different example

## Files 


## Functions 


## Data application
Code related to fitting model to Medicare database is found in data_application

`0_data_prep.R` is script to load data and prepare it for analysis. This is run each time a new input dataset is prepared. This also prepares bootstrap samples
`0_tune_CausalGPS.R` is script that tunes the XGBoost and caliper hyperparameters for the CausalGPS package, this only needs to be run once for each input dataset to find optimal parameters
`1_model_fit.R` is script to fit all ERC estimators to the data and generate the ERC for the MEdicare dataset
`1_model_fit_boot.R` is script to fit all ERC estimators to the bootstrap data and generate the ERC for the Medicare dataset
`2_post_processing.R` is script to generate plots of ERC as well as covariate balance plots with included uncertainty from bootstrap
`/functions/` contains two scripts, the first `/functions/entropy_wt_functions.R` is a more efficient implementation of entropy weighting for a large dataset and `/functions/fit_causalGPS_by_year.R` fits the CausalGPS package by year for Medicare dataset
`bash*.sh` all these scripts are used to launch jobs on the FASSE cluster 
