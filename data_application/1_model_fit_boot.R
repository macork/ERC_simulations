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
library(cobalt)


# Get command arguments (first argument is which input data to use, second is name for model run)
args <- commandArgs(T)
input_flag <- as.character(args[[1]])
model_flag = as.character(args[[2]])
post_only = F

# Define sample number
sample <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/boostrap/boot_data_", sample, ".RDS"))

# Make sure data is in correct form
data <- data.table(data)

# Create output directory for models
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
dir.create(out_dir)

# Now add boostrap component
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/boostrap/")
dir.create(out_dir)

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

# Create logical for running post estimation only 
if (!post_only) {
  # Fit a glm
  linear_fit <- 
    lm(reformulate(c("pm25", strata_var, confounders), response = "log_mort"),
       weight = time_count, data = data)
  saveRDS(linear_fit, file = paste0(out_dir, "linear.RDS"))
    
    # # Fit two gam models for now, one with cubic spline and one with beta spline (penalized)
    gam_fit <- 
      mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
                weight = time_count, data = data)
    saveRDS(gam_fit, file = paste0(out_dir, "gam.RDS"))
    
    # Now fit the entropy based weights in your model 
    source(paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/functions/entropy_wt_functions.R"))
    
    # Fitting entropy weights on zip code level variables
    data_ent <- data %>% select(all_of(c("zip", "pm25", confounders)))
    
    # Get rid of repeats since eliminating strata to fit the models
    data_ent <- unique(data_ent)
    
    # Create model matrx
    data_ent_matrix <- model.matrix(reformulate(c("-1", confounders)),
                                    data = data_ent)
    #c_mat <- scale(data_ent_matrix, scale = apply(data_ent_matrix, 2, function(x) 0.5*diff(range(x))))
    #e <- ebal(data_ent$pm25_ensemble, c_mat)
    
    
    # First center non-binary variables 
    # c_mat <- scale(data_ent_matrix,
    #                center = apply(c_mat, 2, function(x) {
    #                  ifelse(length(unique(x)) > 2, mean(x), 0)
    #                }),
    #                scale = FALSE)
    

    #c_mat <- scale(data_ent_matrix, scale = apply(data_ent_matrix, 2, function(x) 0.5*diff(range(x))))
    # Now scale to be between -1 and 1
    c_mat <- scale(data_ent_matrix,
                   center = T,
                   scale = apply(data_ent_matrix, 2, function(x) ifelse(max(abs(x)) == 0, 1, max(abs(x)))))
    
    #c_mat <- scale(data_ent_matrix, scale = apply(data_ent_matrix, 2, function(x) 0.5*diff(range(x))))
    e <- ebal(data_ent$pm25, c_mat)
    ent_weights <- e$weights
    ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
    
    # # Run through one more iteration to try to get rid of extreme weights
    e <- ebal(data_ent$pm25, c_mat, base_weights = ent_weights)
    ent_weights <- e$weights
    
    # Now save entropy weights and covariate weights
    
    #cov.wt(cbind(data_ent$pm25, c_mat), ent_weights, TRUE)$cor[, 1]
    #cov.wt(cbind(data_ent$pm25, data_ent_matrix), ent_weights, TRUE)$cor[, 1]
    # ent_dataset <- data.table(cbind(pm25_ensemble = data_ent$pm25_ensemble, data_ent_matrix, ent_weights))
    
    # Add the entropy weights to this dataset
    data_ent$ent_weight <- ent_weights
    
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
    #saveRDS(ent_weights, file = paste0(out_dir, "entropy_weights_raw.RDS"))
    
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
    
    
    # Fit two gam models for now, one with cubic spline and one with beta spline (penalized)
    gam_ent <-
      mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
                weight = ent_pop_weight, data = data)
    
    saveRDS(gam_ent, file = paste0(out_dir, "/gam_ent.RDS"))
    
    
    # gam_ent_bs <- 
    #   mgcv::bam(log_mort ~ s(pm25_ensemble, bs = 'bs', k = 20) + medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +as.factor(year) + as.factor(region) + 
    #               as.factor(sex) + as.factor(race) + as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year),
    #             data = data, weights = ent_weights)
    # 
    # saveRDS(gam_ent_bs, file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/gam_ent_bs.RDS"))
    
    # # USe MGCV bam instead of gam 
    # change_ent <-
    #   chngptm(log_mort ~ medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
    #             as.factor(year) + as.factor(region) + as.factor(sex) + as.factor(race) + as.factor(dual) + 
    #             as.factor(entry_age_break) + as.factor(followup_year), ~ pm25_ensemble,
    #           data = data, family = "gaussian", 
    #           type = "segmented", var.type = "default", weights = ent_weights)
    # 
    # saveRDS(change_ent, file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/change_ent.RDS"))
  
  
  # Remember this for fitting only causalGPS package
  # data_ent <-
  #   data[, c("medhouseholdincome", "medianhousevalue", "poverty", "education", "pct_owner_occ",
  #            "year", "region", "sex", "race", "dual", "entry_age_break", "followup_year")]
  
    # Keep only confounders (fiting to zip code, year data)
    data_confounders <- data.frame(data_ent %>% select(all_of(confounders))) %>% select(-year)
    
    # data_ent[, year := as.numeric(as.character(year))]
    source(paste0(proj_dir, "/ERC_simulation/Simulation_studies/functions/data_application_functions.R"))
    
    # Now fit CausalGPS model
    grid_min_results <- 
      fit_causalGPS_by_year(nrounds = 100, 
                            max_depth = 5,
                            eta = 0.3, 
                            delta = 2, 
                            return_corr_only = F, 
                            nthread = 10, trim = T,
                            data_ent = data_ent,
                            data_confounders = data_confounders)
    
    pseudo_pop <- grid_min_results$pseudo_pop
    
    saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
    
  # Now fit semi-parametric here after generating causalGPS 
  #pseudo_pop <- readRDS(file = paste0(out_dir, "pseudo_pop.RDS"))
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
  
  # Multiply weights and use 
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

relative_rate_data <- 
  data_prediction %>% 
  mutate(linear_model = linear_model - data_prediction$linear_model[1]) %>%
  mutate(linear_ent = linear_ent - data_prediction$linear_ent[1]) %>%
  mutate(gam_model = gam_model - data_prediction$gam_model[1]) %>%
  mutate(gam_ent = gam_ent - data_prediction$gam_ent[1]) %>%
  mutate(causal_model = causal_model - data_prediction$causal_model[1]) %>%
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25) %>%
  mutate(prediction = exp(prediction))

saveRDS(relative_rate_data, file = paste0(out_dir, "relative_rate.RDS"))

relative_rate_plot <-
  relative_rate_data %>%
  ggplot(aes(x = pm25, y = exp(prediction), color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Relative mortality rate")

ggsave(relative_rate_plot, file = paste0(out_dir, "relative_rate_plot.pdf"))

message("done")

# # Plot the contrasts
# contrasts <- c(8, 9, 10, 12)
# 
# contrast_prediction <- 
#   rbindlist(lapply(contrasts, function(pot_exp) {
#     
#     potential_data <- mutate(data, pm25_ensemble = pot_exp)
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
#       dplyr::select(pm25_ensemble, linear_model, linear_ent, gam_model, gam_ent, causal_model) %>% 
#       summarize_all(mean) %>% 
#       data.table()
#     return(potential_outcome)
#   }))
# 
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
# # Bind together and save
# contrast_data <- 
#   rbind(contrast8_12, contrast9_12, contrast10_12) %>% 
#   cbind(contrast = c(8, 9, 10)) %>% 
#   data.table() 
# 
# saveRDS(contrast_data, file = paste0(out_dir, "contrast_data.RDS"))




