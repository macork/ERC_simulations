rm(list = ls())

# Load required packages
library(data.table)
library(CausalGPS)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(WeightIt)
library(mgcv)
library(Rcpp)
library(RcppEigen)
library(cobalt)


# Get command arguments (first argument is which input data to use, second is name for model run)
args <- commandArgs(T)
input_flag <- as.character(args[[1]])
model_flag = as.character(args[[2]])
post_only = F

# Define sample number
sample <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load in data and set out directory
proj_dir <- "~/nsaph_projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/boostrap/boot_data_", sample, ".RDS"))

# Make sure data is in correct form
data <- data.table(data)

# Create output directory for models
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
dir.create(out_dir)

# Now add boostrap component
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/boostrap/")
dir.create(out_dir)

# out directory is for each bootstrap sample
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/boostrap/", sample, "/")
dir.create(out_dir)

# Assign variables as strata or confounders
strata_var <- c("female", "race", "dual", "entry_age_break", "followup_year")
confounders <- c("smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "mean_bmi",
                 "poverty", "education", "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax",
                 "winter_rmax", "year", "region")


# Grab min and max for making curve for bootstrap (to make sure same as other dataset)
min_max_curve <- function(input_flag) {
  entire_data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/",
                         input_flag, "/input_data.RDS"))
  # find min and max
  min_pm <- min(entire_data$pm25)
  max_pm <- max(entire_data$pm25)
  return(c(min_pm, max_pm))
}
min_max <- min_max_curve(input_flag)
min_pm <- min_max[1]
max_pm <- min_max[2]

# Run through model fit
if (!post_only) {
  # Fit linear model
  linear_fit <- 
    lm(reformulate(c("pm25", strata_var, confounders), response = "log_mort"),
       weight = time_count, data = data)
  saveRDS(linear_fit, file = paste0(out_dir, "linear.RDS"))
    
  # Fit gam
  gam_fit <- 
    mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
              weight = time_count, data = data)
  saveRDS(gam_fit, file = paste0(out_dir, "gam.RDS"))
    
  # Now fit the entropy based weights in your model 
  #source(paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/functions/entropy_wt_functions.R"))
    
  # Fitting entropy weights on zip code level variables
  data_ent <- data %>% select(all_of(c("zip", "pm25", confounders)))
    
  # Get rid of repeats since eliminating strata to fit the models
  data_ent <- unique(data_ent)
  
  # Try weightit function to calculate entropy balancing weights with second moment
  ent_weights <- 
    weightit(reformulate(confounders, response = "pm25"), data = data_ent,
             method = "ebal", stabilize = TRUE, moments = 2)
  
  # Extract weights
  ent_weights <- ent_weights$weights
    
  # Create model matrx
  data_ent_matrix <- model.matrix(reformulate(c("-1", confounders)),
                                    data = data_ent)
    
  # Now scale to be between -1 and 1
  c_mat <- scale(data_ent_matrix,
                 center = T,
                 scale = apply(data_ent_matrix, 2, function(x) ifelse(max(abs(x)) == 0, 1, max(abs(x)))))
    
  # e <- ebal(data_ent$pm25, c_mat)
  # ent_weights <- e$weights
  
  # Truncate to 99.5% percentile
  ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
    
  # # Run through one more iteration to try to get rid of extreme weights
  # e <- ebal(data_ent$pm25, c_mat, base_weights = ent_weights)
  # ent_weights <- e$weights
    
  # Add the entropy weights to this dataset
  data_ent$ent_weight <- ent_weights
    
  # Make sure balanced
  ent_balance <- 
    bal.tab(reformulate(confounders, response = "pm25"),
            data = data_ent,
            weights = data_ent[["ent_weight"]],
            method = "weighting",
            stats = c("cor"),
            un = T,
            continuous = "std",
            s.d.denom = "weighted",
            abs = T,
            thresholds = c(cor = .1), poly = 1)
    
    # Save entropy weight output
    saveRDS(data_ent, file = paste0(out_dir, "entropy_weights.RDS"))
    
    # Now join entropy weights to original dataset (given you matched on zip and year only)
    data_ent_join <- data_ent[, c("zip", "year", "ent_weight")]
    
    data <- 
      data %>% 
      left_join(data_ent_join, by = c("zip", "year")) %>%
      mutate(ent_pop_weight = time_count * ent_weight) # Multiply survey weights
    
    # Now fit the entropy based weights 
    linear_ent <- 
      lm(reformulate(c("pm25", strata_var, confounders), response = "log_mort"),
         weight = ent_pop_weight, data = data)
    
    saveRDS(linear_ent, file = paste0(out_dir, "linear_ent.RDS"))
    
    
    # Fit GAM with entropy weights
    gam_ent <-
      mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
                weight = ent_pop_weight, data = data)
    
    saveRDS(gam_ent, file = paste0(out_dir, "/gam_ent.RDS"))
    
    # Fit CausalGPS 
    # Keep only confounders (fiting to zip code, year data)
    data_confounders <- data.frame(data_ent %>% select(all_of(confounders))) %>% select(-year)
    
    # Source function to fit causalGPS
    source(paste0(proj_dir, "/ERC_simulation/Simulation_studies/data_application/functions/fit_causalGPS_by_year.R"))
    
    # Now fit CausalGPS model
    grid_min_results <-
      fit_causalGPS_by_year(nrounds = 100,
                            max_depth = 5,
                            eta = 0.3,
                            delta = 2,
                            return_corr_only = F,
                            nthread = 10, 
                            trim = T,
                            data_ent = data_ent,
                            confounders = confounders)
    
    pseudo_pop <- grid_min_results$pseudo_pop
    saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
    
  
    # Join causal weights
    pseudo_pop_frame <- data.table(pseudo_pop)

  # Create causalGPS join to add weights from design stage 
  causalgps_join <-
    pseudo_pop_frame %>%
    select(zip, year, causal_weight = counter_weight)

  data_causalGPS <-
    data %>%
    left_join(causalgps_join, by = c("zip", "year")) %>%
    mutate(causal_pop_weight = time_count * causal_weight) %>% # Multiply survey weights
    data.table()

  # Multiply weights and use to fit causal model
  causal_gps <-
    mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders, "year"), response = "log_mort"),
              weight = causal_pop_weight, data = data_causalGPS)
  saveRDS(causal_gps, file = paste0(out_dir, "causal_fit.RDS"))
  
  message("Done with fitting models")
}

# Load model fits if running post only
linear_fit <- readRDS(file = paste0(out_dir, "linear.RDS"))
gam_fit <- readRDS(file = paste0(out_dir, "gam.RDS"))
linear_ent <- readRDS(file = paste0(out_dir, "linear_ent.RDS"))
gam_ent <- readRDS(file = paste0(out_dir, "gam_ent.RDS"))
causal_gps <- readRDS(file = paste0(out_dir, "causal_fit.RDS"))

# Now get potential outcomes
data_prediction <- 
  rbindlist(lapply(seq(min_pm, max_pm, length.out = 100), function(pot_exp) {
    
    # Get potential data if all had same potential exposure
    potential_data <- 
      select(data, all_of(c(strata_var, confounders))) %>% 
      mutate(pm25 = pot_exp)
    
    # Fit each model and take mean for ERC
    potential_outcome <- 
      potential_data %>% 
      mutate(
        linear_model = predict(linear_fit, potential_data, type = "response"),
        gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
        linear_ent = predict(linear_ent, potential_data, type = "response"),
        gam_ent = predict(gam_ent, newdata = potential_data, type = "response"),
        causal_model =  predict(causal_gps, newdata = potential_data, type = "response")
      ) %>% 
      dplyr::select(pm25, linear_model, linear_ent, gam_model, gam_ent, causal_model) %>% 
      summarize_all(mean) %>% 
      data.table()
    return(potential_outcome)
  }))

saveRDS(data_prediction, file = paste0(out_dir, "data_prediction.RDS"))

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

mortality_data <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25) %>% 
  mutate(prediction = exp(prediction))
saveRDS(mortality_data, file = paste0(out_dir, "mortality_data.RDS"))

mortality_plot <- 
  mortality_data %>%
  ggplot(aes(x = pm25, y = prediction, color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Mortality rate")
ggsave(mortality_plot, file = paste0(out_dir, "mortality_plot.pdf"))

message("done")