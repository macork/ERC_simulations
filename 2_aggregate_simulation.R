# Script to run after finished running through simulation
# Aggregates all fo the simulations and stores them in one RDS file

library(tidyverse)
library(data.table)

exp_relationship = "threshold"
adjust_confounder = T
time_stamp = "first_draft2"
replicates = 100

if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

results_dir <- paste0(repo_dir, "/results/", exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")



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

saveRDS(list(metrics = metrics, predictions = predictions, correlation = correlation),
        file = paste0(results_dir, "aggregated_results.RDS"))

# Now run figures code