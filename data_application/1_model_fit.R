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
input_flag <- as.character(args[[1]]) # Usually kevin_trim_90 for this case
model_flag = as.character(args[[2]]) # Where model will be saved

# Flags for whether to fit only CausalGPS package or only post estimation (models already fit)
only_causalgps = F
only_post = T

# load in data
proj_dir <- "~/nsaph_projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/input_data.RDS"))

# Create output directory for models
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
dir.create(out_dir)
model_dir <- out_dir # old naming convention

# Assign variables as strata or confounders
strata_var <- c("female", "race", "dual", "entry_age_break", "followup_year")
confounders <- c("smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "mean_bmi",
                 "poverty", "education", "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax",
                 "winter_rmax", "year", "region")

# Fit a glm
if (!only_post) {
  if (!only_causalgps) {
    linear_fit <- 
      lm(reformulate(c("pm25", strata_var, confounders), response = "log_mort"),
         weight = time_count, data = data)
    
    saveRDS(linear_fit, file = paste0(out_dir, "linear.RDS"))
    
    # Fit GAM model
    gam_fit <- 
      mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
                weight = time_count, data = data)
    
    saveRDS(gam_fit, file = paste0(out_dir, "gam.RDS"))
    
    
    # Change point model did not work with current dataset
    
    # Now fit entropy weighting 
    # Using custom entropy weighting for first moment
    # Switched to WeightIt package for second moment 
    source(paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/functions/entropy_wt_functions.R"))
    
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
    
    # Create data entropy matrix
    data_ent_matrix <- model.matrix(reformulate(c("-1", confounders)),
                                    data = data_ent)
    
    # Now scale to be between -1 and 1 for model to fit
    c_mat <- scale(data_ent_matrix,
                   center = T,
                   scale = apply(data_ent_matrix, 2, function(x) ifelse(max(abs(x)) == 0, 1, max(abs(x)))))
    
    # # Run entropy weighting algorithm
    # e <- ebal(data_ent$pm25, c_mat)
    # ent_weights <- e$weights
    
    # Truncate 99.5% to lessen extreme weights
    ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
    
    # Run through one more iteration to try to get rid of extreme weights
    # e <- ebal(data_ent$pm25, c_mat, base_weights = ent_weights)
    # ent_weights <- e$weights
    
    # Make sure balanced
    cov.wt(cbind(data_ent$pm25, c_mat), ent_weights, TRUE)$cor[, 1]
    cov.wt(cbind(data_ent$pm25, data_ent_matrix), ent_weights, TRUE)$cor[, 1]
    # ent_dataset <- data.table(cbind(pm25_ensemble = data_ent$pm25_ensemble, data_ent_matrix, ent_weights))
    
    # Add the entropy weights to this dataset
    data_ent$ent_weight <- ent_weights
    
    # Save entropy weight output
    saveRDS(data_ent, file = paste0(out_dir, "entropy_weights.RDS"))
    
    # Now join entropy weights to original dataset (given you matched on zip and year only)
    data_ent_join <- data_ent[, c("zip", "year", "ent_weight")]
    
    # Join and multiply the weights
    data <- 
      data %>% 
      left_join(data_ent_join, by = c("zip", "year")) %>%
      mutate(ent_pop_weight = time_count * ent_weight) # Multiply survey weights
    
    # Now fit the entropy based weights
    linear_ent <- 
      lm(reformulate(c("pm25", strata_var, confounders), response = "log_mort"),
         weight = ent_pop_weight, data = data)
    saveRDS(linear_ent, file = paste0(out_dir, "linear_ent.RDS"))
    
    
    # Fit gam entropy
    gam_ent <-
      mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
                weight = ent_pop_weight, data = data)
    
    saveRDS(gam_ent, file = paste0(out_dir, "/gam_ent.RDS"))
    
    
    # Fit change point entropy
    # change_ent <-
    #   chngptm(reformulate(c(strata_var, confounders), response = "log_mort"), ~ pm25,
    #           data = data, family = "gaussian", 
    #           type = "segmented", var.type = "default", weights = ent_pop_weight)
    # 
    # saveRDS(change_ent, file = paste0(out_dir, "/change_ent.RDS"))
  } 
  
  # Now fit the CausalGPS package
  # Keep only confounders (fiting to zip code, year data)
  data_confounders <- data.frame(data_ent %>% select(all_of(confounders))) %>% select(-year)
  
  # Source function to fit causalGPS by year
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
                          data_confounders = data_confounders)
  
  # Extract pseudo-pop with matching weight and save
  pseudo_pop <- grid_min_results$pseudo_pop
  saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
  
  # Add matching weight to full data set
  psuedo_pop_frame <- data.table(pseudo_pop)

  causalgps_join <- 
    data_ent %>% 
    select(zip, year) %>%
    mutate(causal_weight = psuedo_pop_frame$counter_weight)
  
  # Multiply weights
  data <- 
    data %>% 
    left_join(causalgps_join, by = c("zip", "year")) %>%
    mutate(causal_pop_weight = time_count * causal_weight) # Multiply survey weights
  
  # Now fit causal model with strata specific variables in the outcome model
  causal_gps <- 
    mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
              weight = causal_pop_weight, data = data)
  
  saveRDS(causal_gps, file = paste0(out_dir, "causal_fit.RDS"))
  
  message("Done with fitting models")
}

# load the models if only doing post processing 
linear_fit <- readRDS(file = paste0(out_dir, "linear.RDS"))
gam_fit <- readRDS(file = paste0(out_dir, "gam.RDS"))
linear_ent <- readRDS(file = paste0(out_dir, "linear_ent.RDS"))
gam_ent <- readRDS(file = paste0(out_dir, "gam_ent.RDS"))

# Load causal fit
pseudo_pop <- readRDS(file = paste0(out_dir, "pseudo_pop.RDS"))
causal_gps <- readRDS(file = paste0(out_dir, "causal_fit.RDS"))

# Plot curve from in to max
min_pm <- min(data$pm25)
max_pm <- max(data$pm25)

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
        #change_model = predict(change_model, newdata = potential_data, type = "response"),
        linear_ent = predict(linear_ent, potential_data, type = "response"),
        gam_ent = predict(gam_ent, newdata = potential_data, type = "response"),
        #change_ent = predict(change_ent, newdata = potential_data, type = "response"),
        causal_model =  predict(causal_gps, newdata = potential_data, type = "response")
      ) %>% 
      dplyr::select(pm25, linear_model, linear_ent, gam_model, gam_ent, causal_model) %>% 
      summarize_all(mean) %>% 
      data.table()
    return(potential_outcome)
  }))

saveRDS(data_prediction, file = paste0(out_dir, "data_prediction.RDS"))
data_prediction <- readRDS(file = paste0(out_dir, "data_prediction.RDS"))

# Specify model types
model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

# Create mortality estimates by model type
mortality_data <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25) %>% 
  mutate(prediction = exp(prediction))
saveRDS(mortality_data, file = paste0(out_dir, "mortality_data.RDS"))

# Create plot of ERC by model
mortality_plot <- 
  mortality_data %>%
  ggplot(aes(x = pm25, y = prediction, color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Mortality rate")
ggsave(mortality_plot, file = paste0(out_dir, "mortality_plot.pdf"))

message("done")

# Plot the contrasts (not undertaken for now)
# contrasts <- c(8, 9, 10, 12)
# 
# contrast_prediction <- 
#   rbindlist(lapply(contrasts, function(pot_exp) {
#     
#     potential_data <- mutate(data, pm25 = pot_exp)
#     
#     # Fit each model and take mean for ERC
#     potential_outcome <- 
#       potential_data %>% 
#       mutate(
#         linear_model = predict(linear_fit, potential_data, type = "response"),
#         gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
#         #change_model = predict(change_model, newdata = potential_data, type = "response"),
#         linear_ent = predict(linear_ent, potential_data, type = "response"),
#         gam_ent = predict(gam_ent, newdata = potential_data, type = "response"),
#         #change_ent = predict(change_ent, newdata = potential_data, type = "response"),
#         causal_model =  predict(causal_gps, newdata = rename(potential_data, w = pm25_ensemble), type = "response")
#       ) %>% 
#       dplyr::select(pm25, linear_model, linear_ent, gam_model, gam_ent, causal_model) %>% 
#       summarize_all(mean) %>% 
#       data.table()
#     return(potential_outcome)
#   }))


# contrast8_12 <- 
#   1 - exp(contrast_prediction[pm25_ensemble == 8] - contrast_prediction[pm25_ensemble == 12]) %>%
#   select(!!model_types)
# 
# contrast9_12 <- 
#   1 - exp(contrast_prediction[pm25_ensemble == 9] - contrast_prediction[pm25_ensemble == 12]) %>%
#   select(!!model_types)
# 
# contrast10_12 <- 
#   1 - exp(contrast_prediction[pm25_ensemble == 10] - contrast_prediction[pm25_ensemble == 12]) %>%
#   select(!!model_types)
# 
# contrast_12 <- 
#   contrast_prediction[pm25_ensemble == 12]  %>%
#   select(!!model_types)
# 
# # Bind together and save
# contrast_data <- 
#   rbind(contrast8_12, contrast9_12, contrast10_12, contrast_12) %>% 
#   cbind(contrast = c(8, 9, 10, 12)) %>% 
#   data.table() 
# 
# saveRDS(contrast_data, file = paste0(out_dir, "contrast_data.RDS"))