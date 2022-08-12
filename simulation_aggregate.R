# Script to run after finished running through simulation

library(tidyverse)
library(data.table)

exp_relationship = "sublinear"
adjust_confounder = T
time_stamp = "sublinear_interaction_complex_small"
replicates = 100

results_dir <- paste0("/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/results/",
                      exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")



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