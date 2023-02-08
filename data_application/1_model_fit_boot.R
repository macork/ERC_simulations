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


# Get command arguments (first argument is which input data to use, second is name for model run)
args <- commandArgs(T)
input_flag <- as.character(args[[1]])
model_flag = as.character(args[[2]])
only_causalgps = F

sample <- 1
# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/boostrap/boot_data_", sample, ".RDS"))

# Make sure data is in correct form
data <- data.table(data)
data[, year := factor(year)]
data[, followup_year := factor(followup_year)]

# Create output directory for models
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
dir.create(out_dir)
# Now add boostrap component
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/boostrap/")
dir.create(out_dir)

out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/boostrap/", sample, "/")
dir.create(out_dir)

# Fit a glm
  linear_fit <-
    lm(log_mort ~ pm25_ensemble + medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +as.factor(year) + as.factor(region) +
         as.factor(sex) + as.factor(race) + as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year),
       data = data)
  saveRDS(linear_fit, file = paste0(out_dir, "linear.RDS"))
  
  # # Fit two gam models for now, one with cubic spline and one with beta spline (penalized)
  gam_fit <-
    mgcv::bam(log_mort ~ s(pm25_ensemble, bs = 'cr', k = 4) + medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +as.factor(year) + as.factor(region) +
                as.factor(sex) + as.factor(race) + as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year),
              data = data)
  saveRDS(gam_fit, file = paste0(out_dir, "gam.rds"))
  
  # Now fit the entropy based weights in your model 
  source(paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/functions/entropy_wt_functions.R"))
  
  data_ent <-
    data[, c("medhouseholdincome", "medianhousevalue", "poverty", "education",
             "pct_owner_occ", "year", "region", "sex", "race", "dual", "entry_age_break")]
  
  # Create model matrx
  data_ent_matrix <- model.matrix( ~ -1 + medhouseholdincome + medianhousevalue + poverty + education +
                                     pct_owner_occ + year + region + sex + race + dual + entry_age_break, 
                                   data = data_ent)
  
  #c_mat <- scale(data_ent_matrix, scale = apply(data_ent_matrix, 2, function(x) 0.5*diff(range(x))))
  #e <- ebal(data_ent$pm25_ensemble, c_mat)
  
  
  # First center non-binary variables 
  # c_mat <- scale(data_ent_matrix,
  #                center = apply(c_mat, 2, function(x) {
  #                  ifelse(length(unique(x)) > 2, mean(x), 0)
  #                }),
  #                scale = FALSE)
  
  
  # Now scale to be between -1 and 1
  c_mat <- scale(data_ent_matrix,
                 center = T,
                 scale = apply(data_ent_matrix, 2, function(x) ifelse(max(abs(x)) == 0, 1, max(abs(x)))))
  
  #c_mat <- scale(data_ent_matrix, scale = apply(data_ent_matrix, 2, function(x) 0.5*diff(range(x))))
  e <- ebal(data$pm25_ensemble, c_mat)
  ent_weights <- e$weights
  ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
  
  # # Run through one more iteration to get rid of extreme weights
  e <- ebal(data$pm25_ensemble, c_mat, base_weights = ent_weights)
  ent_weights <- e$weights
  ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
  
  # Now save entropy weights and covariate weights
  cov.wt(cbind(data$pm25_ensemble, c_mat), ent_weights, TRUE)$cor[, 1]
  cov.wt(cbind(data$pm25_ensemble, data_ent_matrix), ent_weights, TRUE)$cor[, 1]
  ent_dataset <- data.table(cbind(pm25_ensemble = data$pm25_ensemble, data_ent_matrix, ent_weights))
  
  # Save entropy weight output
  saveRDS(ent_dataset, file = paste0(out_dir, "entropy_weights.RDS"))
  #saveRDS(ent_weights, file = paste0(out_dir, "entropy_weights_raw.RDS"))
  
  # Now fit the entropy based weights 
  linear_ent <-
    lm(log_mort ~ pm25_ensemble + medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
         as.factor(year) + as.factor(region) +
         as.factor(sex) + as.factor(race) + as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year),
       data = data, weights = ent_weights)
  
  saveRDS(linear_ent, file = paste0(out_dir, "linear_ent.RDS"))
  
  
  # Fit two gam models for now, one with cubic spline and one with beta spline (penalized)
  gam_ent <-
    mgcv::bam(log_mort ~ s(pm25_ensemble, bs = 'cr', k = 4) + medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +as.factor(year) + as.factor(region) +
                as.factor(sex) + as.factor(race) + as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year),
              data = data, weights = ent_weights)
  
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

data_ent <-
  data[, c("medhouseholdincome", "medianhousevalue", "year", "poverty", "education",
           "pct_owner_occ", "region", "sex", "race", "dual", "entry_age_break")]

# data_ent[, year := as.numeric(as.character(year))]

# Now fit CausalGPS model
pseudo_pop <- 
  generate_pseudo_pop(data$log_mort,
                      data$pm25_ensemble,
                      data.frame(data_ent),
                      ci_appr = "matching",
                      pred_model = "sl",
                      gps_model = "parametric",
                      use_cov_transform = TRUE,
                      transformers = list("pow2", "pow3", "abs", "scale"),
                      trim_quantiles = c(0, 1),
                      optimized_compile = TRUE,
                      sl_lib = c("m_xgboost"),
                      params = list(xgb_nrounds = c(50),
                                    xgb_eta = c(.3),
                                    max_depth = c(6),
                                    "xgb_min_child_weight" = 1),
                      covar_bl_method = "absolute",
                      covar_bl_trs = 1,
                      covar_bl_trs_type = "mean",
                      max_attempt = 1,
                      matching_fun = "matching_l1",
                      delta_n = 0.16,
                      scale = 1,
                      nthread = 40)
# If passed threshold keep, if not try other values
if (pseudo_pop$adjusted_corr_results$mean_absolute_corr < 0.1) {
  saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
} else {
  # Now run again with defaults if default is not satisfied
  pseudo_pop2 <- 
    generate_pseudo_pop(data$log_mort,
                        data$pm25_ensemble,
                        data.frame(data_ent),
                        ci_appr = "matching",
                        pred_model = "sl",
                        gps_model = "parametric",
                        use_cov_transform = TRUE,
                        transformers = list("pow2", "pow3", "abs", "scale"),
                        trim_quantiles = c(0, 1),
                        optimized_compile = TRUE,
                        sl_lib = c("m_xgboost"),
                        params = list(xgb_nrounds = c(10, 20, 50, 200),
                                      xgb_eta = c(0.025, 0.05, 0.1, 0.2, 0.3),
                                      max_depth = c(4, 5, 6)),
                        covar_bl_method = "absolute",
                        covar_bl_trs = 0.1,
                        covar_bl_trs_type = "mean",
                        max_attempt = 40,
                        matching_fun = "matching_l1",
                        delta_n = 0.16,
                        scale = 1,
                        nthread = 40)
  
  # If passes then go ahead or load other fit
  if (pseudo_pop2$passed_covar_test) {
    saveRDS(pseudo_pop2, file = paste0(out_dir, "pseudo_pop.RDS"))
  } else {
    saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
  }
}


# Now fit semi-parametric here after generating causalGPS 
pseudo_pop <- readRDS(file = paste0(out_dir, "pseudo_pop.RDS"))
psuedo_pop_frame <- data.table(pseudo_pop$pseudo_pop)

causal_gps <- 
  mgcv::bam(formula = Y ~ s(w, bs = 'cr', k = 4),
            family = "gaussian",
            data = psuedo_pop_frame,
            weights = counter_weight)

saveRDS(causal_gps, file = paste0(out_dir, "causal_fit.RDS"))

message("Done with fitting models")

# Now aggregate data in post processing step

# find min and max
min_pm <- min(data$pm25_ensemble)
max_pm <- max(data$pm25_ensemble)

data_prediction <- 
  rbindlist(lapply(seq(min_pm, max_pm, length.out = 100), function(pot_exp) {
    
    potential_data <- mutate(data, pm25_ensemble = pot_exp)
    
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
        causal_model =  predict(causal_gps, newdata = rename(potential_data, w = pm25_ensemble), type = "response")
      ) %>% 
      dplyr::select(pm25_ensemble, linear_model, linear_ent, gam_model, gam_ent, causal_model) %>% 
      summarize_all(mean) %>% 
      data.table()
    return(potential_outcome)
  }))

saveRDS(data_prediction, file = paste0(model_dir, "data_prediction.RDS"))

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

mortality_data <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25_ensemble) %>% 
  mutate(prediction = exp(prediction))

saveRDS(mortality_data, file = paste0(out_dir, "mortality_data.RDS"))

mortality_plot <- 
  mortality_data %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
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
  arrange(pm25_ensemble) %>%
  mutate(prediction = exp(prediction))

saveRDS(relative_rate_data, file = paste0(out_dir, "relative_rate.RDS"))

relative_rate_plot <-
  relative_rate_data %>%
  ggplot(aes(x = pm25_ensemble, y = exp(prediction), color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Relative mortality rate")

ggsave(relative_rate_plot, file = paste0(model_dir, "relative_rate_plot.pdf"))

# Plot the contrasts
contrasts <- c(8, 9, 10, 12)

contrast_prediction <- 
  rbindlist(lapply(contrasts, function(pot_exp) {
    
    potential_data <- mutate(data, pm25_ensemble = pot_exp)
    
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
        causal_model =  predict(causal_gps, newdata = rename(potential_data, w = pm25_ensemble), type = "response")
      ) %>% 
      dplyr::select(pm25_ensemble, linear_model, linear_ent, gam_model, gam_ent, causal_model) %>% 
      summarize_all(mean) %>% 
      data.table()
    return(potential_outcome)
  }))


contrast8_12 <- 
  1 - exp(contrast_prediction[pm25_ensemble == 8] - contrast_prediction[pm25_ensemble == 12]) %>%
  select(!!model_types)

contrast9_12 <- 
  1 - exp(contrast_prediction[pm25_ensemble == 9] - contrast_prediction[pm25_ensemble == 12]) %>%
  select(!!model_types)

contrast10_12 <- 
  1 - exp(contrast_prediction[pm25_ensemble == 10] - contrast_prediction[pm25_ensemble == 12]) %>%
  select(!!model_types)

# Bind together and save
contrast_data <- 
  rbind(contrast8_12, contrast9_12, contrast10_12) %>% 
  cbind(contrast = c(8, 9, 10)) %>% 
  data.table() 

saveRDS(contrast_data, file = paste0(out_dir, "contrast_data.RDS"))




