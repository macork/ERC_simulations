# Script to run after finished running through simulation
# Aggregates all of the simulations and stores them in one RDS file

library(tidyverse)
library(data.table)

# set model tag and relationship
exp_relationship = "threshold"
adjust_confounder = T
model_tag = "first_submission4"
replicates = 100

# Grab repo directory (based on either local or on cluster)
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

# Get results directory
results_dir <- paste0(repo_dir, "/results/", exp_relationship, "_", adjust_confounder,
                      "/", model_tag, "/")

# Bind together metrics
metrics <- 
  rbindlist(lapply(1:replicates, function(rep){
    sim_metrics <- try(readRDS(paste0(results_dir, "sim_final_", rep, ".RDS"))[[1]], silent = T)
    if (is.data.frame(sim_metrics)) {
      data.table(sim_metrics)
      return(sim_metrics)
    } else {
      return()
    }
  }))

# Bind together predictions
predictions <- 
  rbindlist(lapply(1:replicates, function(rep){
    sim_pred <- try(readRDS(paste0(results_dir, "sim_final_", rep, ".RDS"))[[2]], silent = T)
    if (is.data.frame(sim_pred)) {
      data.table(sim_pred)
      return(sim_pred)
    } else {
      return()
    }
  }))

# Bind together correlations
correlation <- 
  rbindlist(lapply(1:replicates, function(rep){
    sim_corr <- try(readRDS(paste0(results_dir, "sim_final_", rep, ".RDS"))[[3]], silent = T)
    if (is.data.frame(sim_corr)) {
      data.table(sim_corr)
      return(sim_corr)
    } else {
      return()
    }
  }))

# Save list as aggregated results
saveRDS(list(metrics = metrics, predictions = predictions, correlation = correlation),
        file = paste0(results_dir, "aggregated_results.RDS"))
