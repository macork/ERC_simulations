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
# library("devtools")
#install_github("fasrc/CausalGPS")
library("CausalGPS")


# File paths
# Set directory (based on username on cluster or local computer)
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

# Read in config file and grab arguments
config <- fread(paste0(repo_dir,"config.csv"), header = T)
exp_relationship <- as.character(config[Argument == "exp_relationship", Value])
adjust_confounder <- as.logical(config[Argument == "adjust_confounder", Value])
sample_size <- as.numeric(config[Argument == "sample_size", Value])
out_relationship <- as.character(config[Argument == "out_relationship", Value])
confounder_setting <- as.character(config[Argument == "confounder_setting", Value])

# pass argument that explains run
args <- commandArgs(T)
run_title = as.character(args[[1]])

# Create out directory
out_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder)
if (!dir.exists(out_dir)) dir.create(out_dir)

# Create correct output directory (previously using time, now using whatever character you want)
#out_dir_tag <- format(Sys.time(), "%m%d_%H")
dir.create(paste0(out_dir, "/", run_title))
out_dir <- paste0(out_dir, "/", run_title)

# Write config to output directory for reproducing results
write.csv(config, file = paste0(out_dir, "/config.csv"))

# Source the appropriate functions needed for simulation 
source(paste0(repo_dir, "simulation_functions.R"))

# Define functions used
cov_function <- function(confounders) as.vector(-0.8 + matrix(c(0.1, 0.1, -0.1, 0.2, 0.1, 0.1), nrow = 1) %*% t(confounders))
scale_exposure <- function(x){20 * (x-min(x))/(max(x)-min(x))}

# Define arguments
sim.num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Start the simulations ------------------------------------------------------------------------------------------
set.seed(sim.num)

# Generate the confounders
if (confounder_setting = "simple") {
  cf <- mvrnorm(n = sample_size,
                mu = rep(0, 4),
                Sigma = diag(4))
  cf5 <- sample(c((-2):2), sample_size, replace = T)
  cf6 <- runif(sample_size, min = -3, max = 3)
} else if (confounder_setting == "nonzero") {
  cf <- mvrnorm(n = sample_size,
                mu = rep(2, 4),
                Sigma = diag(4))
  cf5 <- sample(c((-3):2), sample_size, replace = T)
  cf6 <- runif(sample_size, min = -1, max = 4)
} else if (confounder_setting = "correlated") {
  cf <- mvrnorm(n = sample_size,
                mu = rep(0, 4),
                Sigma = diag(x = 0.8, nrow = 4, ncol = 4) + 0.2)
  cf5 <- sample(c((-2):2), sample_size, replace = T)
  cf6 <- runif(sample_size, min = -3, max = 3)
} else if (confounder_setting = "complex"){
  cf <- mvrnorm(n = sample_size,
                mu = rep(2, 4),
                Sigma = diag(x = 0.8, nrow = 4, ncol = 4) + 0.2)
  cf5 <- sample(c((-3):2), sample_size, replace = T)
  cf6 <- runif(sample_size, min = -1, max = 4)
} else {
  stop("Confounder settings are either simple, nonzero, correlated, or complex")
}

# Bing together confounders and name
confounders = cbind(cf, cf5, cf6)
colnames(confounders) = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")

# Loop through all GPS mods
sim_results <- list()

# Loop through each iteration of the GPS model
for (gps_mod in 1:4) {
  # Generate appropriate exposure (make this a function later)
  if (gps_mod == 1) {
    exposure = scale_exposure(scale_exposure(cov_function(confounders)) + rnorm(sample_size, mean = 0, sd = sqrt(10)))
  } else if (gps_mod == 2) {
    exposure = scale_exposure(scale_exposure(cov_function(confounders)) + (sqrt(5)) * rt(sample_size, df = 3))
  } else if (gps_mod == 3) {
    exposure = scale_exposure(scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2) + rnorm(sample_size, mean = 0, sd = sqrt(10)))
  } else if (gps_mod == 4) {
    exposure = scale_exposure(scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2 + 0.5 * (confounders[, "cf1"]) * (confounders[, "cf5"])) + rnorm(sample_size, mean = 0, sd = sqrt(10)))
  }
  
  # Get metrics and predictions from sample
  metrics_predictions <- 
    metrics_from_data(exposure = exposure, confounders = confounders, exposure_relationship = exp_relationship,
                      outcome_relationship = out_relationship, sample_size = sample_size, family = "gaussian", 
                      adjust_confounder = adjust_confounder, causal_gps = T)
  
  metrics <- metrics_predictions$metrics
  predictions <-  metrics_predictions$predictions
  cor_table <- metrics_predictions$cor_table
  
  # Add appropriate columns
  metrics <- metrics %>% mutate(gps_mod = gps_mod, sample_size = sample_size, sim = sim.num)
  predictions <- predictions %>% mutate(gps_mod = gps_mod, sample_size = sample_size, sim = sim.num)
  
  # put method for correlation table
  cor_table <- cor_table %>% mutate(method = "ipw", delta = 0)
  
  # Add on correlation from CausalGPS package implementation
  cor_default <-
    data.table(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
               pre_cor = metrics_predictions[[4]]$original_corr_results$absolute_corr,
               post_cor = metrics_predictions[[4]]$adjusted_corr_results$absolute_corr,
               method = "causal_gps_default",
               delta =  metrics_predictions[[4]]$params$delta_n)
  
  cor_tuned <-
    data.table(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"),
               pre_cor = metrics_predictions[[5]]$original_corr_results$absolute_corr,
               post_cor = metrics_predictions[[5]]$adjusted_corr_results$absolute_corr,
               method = "causal_gps_tuned",
               delta =  metrics_predictions[[5]]$params$delta_n)
  
  cor_table <- rbind(cor_table, cor_default, cor_tuned)
  cor_table <- cor_table %>% mutate(gps_mod = gps_mod, sample_size = sample_size, sim = sim.num)
  
  # Store in list all of the results
  sim_results[[gps_mod]] <- list(data.table(metrics),data.table(predictions), data.table(cor_table))
}

# unpack list of all results
metrics_sim <-rbindlist(lapply(sim_results, function(x) x[[1]]))
predictions_sim <- rbindlist(lapply(sim_results, function(x) x[[2]]))
correlation_sim <- rbindlist(lapply(sim_results, function(x) x[[3]]))

# Save the final sim results in list and write to output directory
sim_final <- list(metrics_sim, predictions_sim, correlation_sim)
saveRDS(sim_final, file = paste0(out_dir, "/sim_final_", sim.num, ".RDS"))
