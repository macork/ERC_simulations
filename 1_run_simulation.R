# Function for running ERC sim on FASSE cluster 
rm(list = ls())

# Load libraries needed 
library(tidyverse)
library(data.table)
library(MASS)
library(parallel)
library(xgboost)
library(SuperLearner)
library(weights)
library(CBPS)
library(caret)
library(WeightIt)
library(chngpt)
library(cobalt)
library("CausalGPS", lib.loc = "/n/home_fasse/mcork/apps/ERC_simulation/R_4.0.5")


# File paths
# Set directory (based on username on cluster or local computer)
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

# Keep adjust confounder set to T
adjust_confounder <- T

# Define simulation number
sim.num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# pass argument for model run
args <- commandArgs(T)
exp_relationship = as.character(args[[1]])
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

# data table to store output of loop
sim_metrics <- data.table()
sim_predictions <- data.table()
sim_cor_table <- data.table()


# loop through all sample sizes included
for (sample_size in c(200, 1000, 10000)) {
  # loop through all types of confounder settings
  for (confounder_setting in c("simple")) {
    # Loop through two different outcome relationships
    for (out_relationship in c("linear", "interaction")) {
    
      # Generate the confounders
      # Over sample sample size to then remove negative values
      if (confounder_setting == "simple") {
        cf <- mvrnorm(n = 2 * sample_size,
                      mu = rep(0, 4),
                      Sigma = diag(4))
        cf5 <- sample(c((-2):2), 2 * sample_size, replace = T)
        cf6 <- runif(2 * sample_size, min = -3, max = 3)
      } else if (confounder_setting == "nonzero") {
        cf <- mvrnorm(n = 2 * sample_size,
                      mu = rep(1, 4),
                      Sigma = diag(4))
        cf5 <- sample(c((-3):2), 2 * sample_size, replace = T)
        cf6 <- runif(2 * sample_size, min = -2, max = 3)
      } else if (confounder_setting == "correlated") {
        cf <- mvrnorm(n = 2 * sample_size,
                      mu = rep(0, 4),
                      Sigma = diag(x = 0.7, nrow = 4, ncol = 4) + 0.3)
        cf5 <- sample(c((-2):2), 2 * sample_size, replace = T)
        cf6 <- runif(2 * sample_size, min = -3, max = 3)
      } else if (confounder_setting == "complex"){
        cf <- mvrnorm(n = 2 * sample_size,
                      mu = rep(2, 4),
                      Sigma = diag(x = 0.7, nrow = 4, ncol = 4) + 0.3)
        cf5 <- sample(c((-3):2), 2 * sample_size, replace = T)
        cf6 <- runif(2 * sample_size, min = -1, max = 4)
      } else {
        stop("Confounder settings are either simple, nonzero, correlated, or complex")
      }
    
      # Bring together confounders and name
      confounders_large = cbind(cf, cf5, cf6)
      colnames(confounders_large) = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")
      
      # Loop through each iteration of the GPS model included 
      for (gps_mod in 1:4) {
        # Generate appropriate exposure
        # Make sure all are positive (and only keep sample size)
        if (gps_mod == 1) {
          x = 9 * cov_function(confounders_large) + 18 + rnorm(sample_size, mean = 0, sd = sqrt(10))
          exposure_df <- 
            cbind(data.frame(confounders_large), exposure = x) %>% 
            filter(exposure > 0) %>% 
            slice_sample(n = sample_size, replace = F)
          exposure = exposure_df$exposure
          confounders = as.matrix(exposure_df[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")])
        } else if (gps_mod == 2) {
          x = 9 * cov_function(confounders_large) + 18 + (sqrt(5)) * rt(sample_size, df = 3)
          exposure_df <- 
            cbind(data.frame(confounders_large), exposure = x) %>% 
            filter(exposure > 0) %>% 
            slice_sample(n = sample_size, replace = F)
          exposure = exposure_df$exposure
          confounders = as.matrix(exposure_df[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")])
        } else if (gps_mod == 3) {
          x = 9 * cov_function(confounders_large) + 15 + 2 * (confounders_large[, "cf3"]) ^ 2 + rnorm(sample_size, mean = 0, sd = sqrt(10))
          exposure_df <- 
            cbind(data.frame(confounders_large), exposure = x) %>%
            filter(exposure > 0) %>% 
            slice_sample(n = sample_size, replace = F)
          exposure = exposure_df$exposure
          confounders = as.matrix(exposure_df[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")])
        } else if (gps_mod == 4) {
          x = 9 * cov_function(confounders_large) + 2 * (confounders_large[, "cf3"]) ^ 2 + 2 * (confounders_large[, "cf1"]) * (confounders_large[, "cf4"]) + 15 + rnorm(sample_size, mean = 0, sd = sqrt(10))
          exposure_df <-
            cbind(data.frame(confounders_large), exposure = x) %>% 
            filter(exposure > 0) %>% 
            slice_sample(n = sample_size, replace = F)
          exposure = exposure_df$exposure
          confounders = as.matrix(exposure_df[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")])
        }
        
        # Get metrics and predictions from sample using function
        metrics_predictions <- 
          metrics_from_data(exposure = exposure, 
                            confounders = confounders, 
                            exposure_relationship = exp_relationship,
                            outcome_relationship = out_relationship, 
                            sample_size = sample_size, 
                            family = "gaussian", 
                            adjust_confounder = adjust_confounder, 
                            causal_gps = T)
        
        metrics <- metrics_predictions$metrics
        predictions <-  metrics_predictions$predictions
        cor_table <- metrics_predictions$cor_table
        
        # Add appropriate columns to metrics and predictions
        metrics <- 
          metrics %>% 
          mutate(gps_mod = gps_mod,
                 sample_size = sample_size,
                 sim = sim.num,
                 confounder_setting = confounder_setting,
                 out_relationship = out_relationship)
        
        predictions <- 
          predictions %>% 
          mutate(gps_mod = gps_mod, 
                 sample_size = sample_size, 
                 sim = sim.num,
                 confounder_setting = confounder_setting,
                 out_relationship = out_relationship)
        
        # put method for correlation table
        cor_table <- cor_table %>% mutate(method = "ipw", delta = 0)
        
        # Add on correlation
        causal_default_pop <- metrics_predictions[[4]]$pseudo_pop
        causal_default_balance <- 
          bal.tab(w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                  data = causal_default_pop,
                  weights = causal_default_pop[["counter_weight"]],
                  method = "weighting",
                  stats = c("cor"),
                  un = T,
                  continuous = "std",
                  s.d.denom = "weighted",
                  abs = T,
                  thresholds = c(cor = .1), poly = 1)
        
        cor_default <-
          data.table(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
                     pre_cor = cor_table$pre_cor,
                     post_cor = causal_default_balance$Balance$Corr.Adj,
                     method = "causal_gps_default",
                     delta =  metrics_predictions[[4]]$params$delta_n)
        
        causal_tuned_pop <- metrics_predictions[[5]]
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
          data.table(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
                     pre_cor = cor_table$pre_cor,
                     post_cor = causal_tuned_balance$Balance$Corr.Adj,
                     method = "causal_gps_tuned",
                     delta =  0)
        
        cor_table <- rbind(cor_table, cor_default, cor_tuned)
        cor_table <- cor_table %>% mutate(gps_mod = gps_mod, sample_size = sample_size, sim = sim.num, 
                                          confounder_setting = confounder_setting, out_relationship = out_relationship)
        
        # Bind together data table with all results
        sim_metrics <- rbind(sim_metrics, data.table(metrics))
        sim_predictions <- rbind(sim_predictions, data.table(predictions))
        sim_cor_table <- rbind(sim_cor_table, data.table(cor_table))
      }
    }
  }
}

# Save the final sim results in list and write to output directory
sim_final <- list(sim_metrics, sim_predictions, sim_cor_table)
saveRDS(sim_final, file = paste0(out_dir, "/sim_final_", sim.num, ".RDS"))
