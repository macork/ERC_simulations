# ERC simulations
Code for running analysis in: Methods for Estimating the Exposure-Response Curve to Inform the New Safety Standards for Fine Particulate Matter

This respository contains code used to conduct simulation studies comparing methods for exposure-response curves. It also contains a data application to the Medicare data set. 

The main code scripts are as follows:
* `0_launch.jobs.sh` is script to launch jobs on the FASSE cluster
* `1_run_simulation.R` is script to launch simulation for one exposure-outcome model, which launches 100 jobs
* `2_aggregate_simulation.R` is script bind together 100 simulation results to calculate predictions and metrics
* `3_produce_figures.R` is script to produce figures to inspect model results
* `4_paper_figures.R` is script to produce figures used in the paper
* `/functions/` contains scripts for running simulation and data application
* `/figures/` saves figures used in publication

## Data application
Code related to fitting model to Medicare database is found in `/data_application/`

* `0_data_prep.R` is script to load data and prepare it for analysis. This is run each time a new input dataset is prepared. This also prepares bootstrap samples
* `0_tune_CausalGPS.R` is script that tunes the XGBoost and caliper hyperparameters for the CausalGPS package, this only needs to be run once for each input dataset to find optimal parameters
* `1_model_fit.R` is script to fit all ERC estimators to the data and generate the ERC for the MEdicare dataset
* `1_model_fit_boot.R` is script to fit all ERC estimators to the bootstrap data and generate the ERC for the Medicare dataset
* `2_post_processing.R` is script to generate plots of ERC as well as covariate balance plots with included uncertainty from bootstrap
* `/functions/` contains two functions used for data application:
  * `/functions/entropy_wt_functions.R` is a more efficient implementation of entropy weighting for a large dataset
  * `/functions/fit_causalGPS_by_year.R` fits the CausalGPS package by year for Medicare dataset
* `bash_*.sh` all these scripts are used to launch jobs on the FASSE cluster 
