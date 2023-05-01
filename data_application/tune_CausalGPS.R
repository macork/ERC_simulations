### File to test covariate balance
rm(list = ls())
# Data application
library(data.table)
library("CausalGPS", lib.loc = "/n/home_fasse/mcork/apps/ERC_simulation/R_4.0.5")
library(tidyverse)
library(chngpt)
library(WeightIt)
library(fastDummies)
library(mgcv)
library(Rcpp)
library(RcppEigen)
library(parallel)
library(cobalt)


# Get command arguments (first argument is which input data to use, 
# second is name for model run)
input_flag <- "kevin_trim_90"
model_flag <- "Causal_Tune_byyear"

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/input_data.RDS"))
#data <- sample_n(data, 1000)
#data[, year := as.factor(year)]

# source function to fit causalGPS
source(paste0(proj_dir, "/ERC_simulation/Simulation_studies/functions/data_application_functions.R"))

# Create output directory for models
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
dir.create(out_dir)
model_dir <- out_dir # old naming convention

# Assign variables as strata or confounders
strata_var <- c("female", "race", "dual", "entry_age_break", "followup_year")
confounders <- c("smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "mean_bmi",
                 "poverty", "education", "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax",
                 "winter_rmax", "region")


data_ent <- data %>% select(all_of(c("zip", "pm25", "year", confounders)))

# Get rid of repeats since eliminating strata to fit the models
data_ent <- unique(data_ent)
data_confounders <- data.frame(data_ent %>% select(all_of(confounders)))

tune_grid <- 
  expand.grid(nrounds = c(50, 100, 150),
              eta = c(0.3, 0.4, 0.45, 0.5),
              max_depth = c(5, 6),
              delta = c(0.2, 0.5, 1, 2))
#tune_grid <- sample_n(tune_grid, 4)

tune_list <- split(tune_grid, seq(nrow(tune_grid)))

grid_search_corr <- 
  rbindlist(mclapply(tune_list, mc.cores = 5, function(x){
    results <- 
      fit_causalGPS_by_year(nrounds = x$nrounds,
                            max_depth = x$max_depth,
                            eta = x$eta,
                            delta = x$delta,
                            return_corr_only = T,
                            nthread = 8, 
                            trim = F)
    return(results)
  }))

saveRDS(grid_search_corr, paste0(out_dir, "/grid_search_corr.RDS"))

grid_min <- grid_search_corr %>% filter(corr == min(corr))
grid_min_results <- 
  fit_causalGPS_by_year(nrounds = grid_min$nrounds, max_depth = grid_min$depth,
                        eta = grid_min$eta, delta = grid_min$delta, return_corr_only = F, 
                        nthread = 10, trim = F)

pseudo_pop <- grid_min_results$pseudo_pop

saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
message("Done")

# load results
grid_search_corr <-  readRDS(paste0(out_dir, "/grid_search_corr.RDS"))
pseudo_pop <- readRDS(file = paste0(out_dir, "pseudo_pop.RDS"))
pseudo_pop_frame <- data.table(pseudo_pop)

# Create causalGPS join (make sure no trimming going forward)
causalgps_join <- 
  pseudo_pop_frame %>% 
  select(zip = Y, year, causal_weight = counter_weight)

data_causalGPS <- 
  data %>% 
  left_join(causalgps_join, by = c("zip", "year")) %>%
  mutate(causal_pop_weight = time_count * causal_weight) %>% # Multiply survey weights
  data.table()

non_causal <- 
  mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders, "year"), response = "log_mort"), data = data_causalGPS)

# Counter only
causal_gps1 <- 
  mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders, "year"), response = "log_mort"),
            weight = causal_weight, data = data_causalGPS)

# Multiply weights
causal_gps2 <- 
  mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders, "year"), response = "log_mort"),
            weight = causal_pop_weight, data = data_causalGPS)

# Counter and strata only
causal_gps3 <- 
  mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var), response = "log_mort"),
            weight = causal_weight, data = data_causalGPS)

# Multiply weights and strata only
causal_gps4 <- 
  mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var), response = "log_mort"),
            weight = causal_pop_weight, data = data_causalGPS)

# Save results
saveRDS(non_causal, file = paste0(out_dir, "non_causal.RDS"))
saveRDS(causal_gps1, file = paste0(out_dir, "causal_gps1.RDS"))
saveRDS(causal_gps2, file = paste0(out_dir, "causal_gps2.RDS"))
saveRDS(causal_gps3, file = paste0(out_dir, "causal_gps3.RDS"))
saveRDS(causal_gps4, file = paste0(out_dir, "causal_gps4.RDS"))

non_causal <- readRDS(file = paste0(out_dir, "non_causal.RDS"))
causal_gps1 <- readRDS(file = paste0(out_dir, "causal_gps1.RDS"))
causal_gps2 <- readRDS(file = paste0(out_dir, "causal_gps2.RDS"))
saveRDS(causal_gps1, file = paste0(out_dir, "causal_gps1.RDS"))
saveRDS(causal_gps2, file = paste0(out_dir, "causal_gps2.RDS"))
saveRDS(causal_gps3, file = paste0(out_dir, "causal_gps3.RDS"))
saveRDS(causal_gps4, file = paste0(out_dir, "causal_gps4.RDS"))
