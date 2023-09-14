# Run sensitivity analysis where data application = T (so taking from data application to run simulation)

# Change library path for running on cluster with updates
.libPaths(new = c("~/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

# Load libraries needed 
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(parallel)
library(xgboost)
library(SuperLearner)
library(WeightIt)
library(chngpt)
library(cobalt)
library(CausalGPS)


# File paths
# Set directory (based on username on cluster or local computer)
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

# Define simulation number
sim.num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# pass argument for model run
args <- commandArgs(T)
exposure_response_relationship = as.character(args[[1]])
run_title = as.character(args[[2]])

# Create out directory
out_dir <- paste0(repo_dir, "results/", exposure_response_relationship, "_TRUE")

# Create correct output directory (previously using time, now using whatever character you want)
dir.create(paste0(out_dir, "/", run_title))
out_dir <- paste0(out_dir, "/", run_title)

# Source the appropriate functions needed for simulation 
source(paste0(repo_dir, "/functions/simulation_functions.R"))

# Load in input data since we are using data application
input_flag <- "kevin_trim_90"
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
out_dir <- paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/")
data_app <- readRDS(file = paste0(out_dir, "/input_data.RDS"))


# Start the simulations ------------------------------------------------------------------------------------------
set.seed(sim.num)

# Store output from for loop
sim_metrics <- tibble()
sim_predictions <- tibble()
sim_cor_table <- tibble()
sim_convergence <- tibble()


# loop through all sample sizes included
for (sample_size in c(200, 1000, 10000)) {
  message(paste("Running with sample size", sample_size))
    # Loop through two different outcome relationships
  for (outcome_interaction in c("T", "F")) {
    if (outcome_interaction) {
      message("Running with interaction in outcome model")
    } else {
      message("Running without interaction in outcome model")
    }
    
    # Create data sample from input data, scale all covariates
    data_app_sample <- 
      data_app %>% 
      sample_n(sample_size, replace = FALSE) %>%
      dplyr::select(exposure = pm25, cf1 = mean_bmi, cf2 = medianhousevalue, 
                    cf3 = smoke_rate, cf4 = education, cf5 = poverty, cf6 = winter_tmmx) %>%
      mutate(across(-exposure, scale))
    
    # Generate synthetic data for simulation study
    sim_data <- sim_data_generate(sample_size = sample_size, 
                                  gps_mod = 1, # not relevant here, placeholder
                                  exposure_response_relationship = exposure_response_relationship, 
                                  outcome_interaction = as.logical(outcome_interaction),
                                  outcome_sd = 10, 
                                  data_application = TRUE,
                                  data_app_sample = data_app_sample)
    
    # Get metrics and predictions from sample
    metrics_predictions <- 
      metrics_from_data(sim_data = sim_data,
                        exposure_response_relationship = exposure_response_relationship, 
                        outcome_interaction = as.logical(outcome_interaction))
    
    metrics <- metrics_predictions$metrics
    predictions <-  metrics_predictions$predictions
    cor_table <- metrics_predictions$cor_table
    convergence_info <- metrics_predictions$convergence_info
    
    # Add appropriate columns to metrics and predictions
    metrics <- 
      metrics %>% 
      mutate(gps_mod = "data_app",
             sample_size = sample_size,
             sim = sim.num,
             exposure_response_relationship = exposure_response_relationship)
    
    predictions <- 
      predictions %>% 
      mutate(gps_mod = "data_app", 
             sample_size = sample_size, 
             sim = sim.num,
             exposure_response_relationship = exposure_response_relationship)
    
    cor_table <- 
      cor_table %>% 
      mutate(gps_mod = "data_app", sample_size = sample_size, sim = sim.num,
             exposure_response_relationship = exposure_response_relationship)
    
    # Add on convergence info for entropy weighting
    convergence_info <- 
      convergence_info %>% 
      mutate(gps_mod = "data_app", sample_size = sample_size, sim = sim.num,
             exposure_response_relationship = exposure_response_relationship)
    
    # Bind together data table with all results
    sim_metrics <- bind_rows(sim_metrics, metrics)
    sim_predictions <- bind_rows(sim_predictions, predictions)
    sim_cor_table <- bind_rows(sim_cor_table, cor_table)
    sim_convergence <- bind_rows(sim_convergence, convergence_info)
  }
}

message("Saving output")

# Save the final sim results in list and write to output directory
sim_final <- list(sim_metrics, sim_predictions, sim_cor_table, sim_convergence)
saveRDS(sim_final, file = paste0(out_dir, "/sim_final_", sim.num, ".RDS"))