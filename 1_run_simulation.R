# Function for running ERC sim on FASSE cluster 
rm(list = ls())

# Load libraries needed 
library(tidyverse)
library(MASS)
library(parallel)
library(xgboost)
library(SuperLearner)
library(weights)
library(CBPS)
library(WeightIt)
library(chngpt)
library(cobalt)
library("CausalGPS")


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
out_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder)
if (!dir.exists(out_dir)) dir.create(out_dir)

# Create correct output directory (previously using time, now using whatever character you want)
dir.create(paste0(out_dir, "/", run_title))
out_dir <- paste0(out_dir, "/", run_title)

# Source the appropriate functions needed for simulation 
source(paste0(repo_dir, "/functions/simulation_functions.R"))

# Define functions used for scaling exposure and gamma function
cov_function <- function(confounders) as.vector(-0.8 + matrix(c(0.1, 0.1, -0.1, 0.2, 0.1, 0.1), nrow = 1) %*% t(confounders))
scale_exposure <- function(x){20 * (x-min(x))/(max(x)-min(x))}

# Start the simulations ------------------------------------------------------------------------------------------
set.seed(sim.num)

# Store output from for loop
sim_metrics <- tibble()
sim_predictions <- tibble()
sim_cor_table <- tibble()


# loop through all sample sizes included
for (sample_size in c(200, 1000, 10000)) {
  # loop through all gps model specifications
  for (gps_mod in 1:4) {
    # Loop through two different outcome relationships
    for (outcome_interaction in c("T", "F")) {
      
      # Generate synthetic data for simulation study
      sim_data <- sim_data_generate(sample_size = 1000, 
                                    gps_mod = gps_mod,
                                    exposure_response_relationship = exposure_response_relationship, 
                                    outcome_interaction = as.logical(outcome_interaction),
                                    outcome_sd = 10, 
                                    data_application = FALSE)
        
        # Get metrics and predictions from sample
        metrics_predictions <- 
          metrics_from_data(sim_data = sim_data,
                            exposure_response_relationship = exp_relationship, 
                            outcome_interaction = as.logical(outcome_interaction))
        
        metrics <- metrics_predictions$metrics
        predictions <-  metrics_predictions$predictions
        cor_table <- metrics_predictions$cor_table
        
        # Add appropriate columns to metrics and predictions
        metrics <- 
          metrics %>% 
          mutate(gps_mod = gps_mod,
                 sample_size = sample_size,
                 sim = sim.num,
                 exposure_response_relationship = exposure_response_relationship)
        
        predictions <- 
          predictions %>% 
          mutate(gps_mod = gps_mod, 
                 sample_size = sample_size, 
                 sim = sim.num,
                 exposure_response_relationship = exposure_response_relationship)
        
        # put method for correlation table
        cor_table <- cor_table %>% mutate(delta = 0)
        
        # Add on correlation from causalGSP
        causal_tuned_pop <- metrics_predictions[[4]]
        causal_tuned_balance <- 
          bal.tab(w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                  data = causal_tuned_pop,
                  weights = causal_tuned_pop[["counter_weight"]],
                  method = "weighting",
                  stats = c("cor"),
                  un = T,
                  continuous = "std",
                  s.d.denom = "weighted",
                  abs = T,
                  thresholds = c(cor = .1), poly = 1)
        
        cor_tuned <-
          data.frame(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
                     pre_cor = cor_table$pre_cor,
                     post_cor = causal_tuned_balance$Balance$Corr.Adj,
                     method = "causal_gps_tuned",
                     delta =  0)
        
        cor_table <- rbind(cor_table, cor_tuned)
        cor_table <- cor_table %>% 
          mutate(gps_mod = gps_mod, sample_size = sample_size, sim = sim.num,
                 exposure_response_relationship = exposure_response_relationship)
        
        # Bind together data table with all results
        sim_metrics <- bind_rows(sim_metrics, metrics)
        sim_predictions <- bind_rows(sim_predictions, predictions)
        sim_cor_table <- bind_rows(sim_cor_table, cor_table)
    }
  }
}

# Save the final sim results in list and write to output directory
sim_final <- list(sim_metrics, sim_predictions, sim_cor_table)
saveRDS(sim_final, file = paste0(out_dir, "/sim_final_", sim.num, ".RDS"))
