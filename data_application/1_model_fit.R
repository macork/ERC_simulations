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
input_flag <- as.character(args[[1]]) # Usually kevin_data for this case
model_flag = as.character(args[[2]])
only_causalgps = F
only_post = T

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/input_data.RDS"))

#data[, year := as.factor(year)]

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
    
    # # Fit two gam models for now, one with cubic spline and one with beta spline (penalized)
    gam_fit <- 
      mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
                weight = time_count, data = data)
    
    saveRDS(gam_fit, file = paste0(out_dir, "gam.RDS"))
    
   
    # 
    # change_model <-
    #   chngptm(log_mort ~ medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
    #             as.factor(year) + as.factor(region) + as.factor(sex) + as.factor(race) + as.factor(dual) + 
    #             as.factor(entry_age_break) + as.factor(followup_year), ~ pm25_ensemble,
    #           data = data, family = "gaussian", 
    #           type = "segmented", var.type = "default")
    # 
    # saveRDS(change_model, file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/change_model.RDS"))
    
    # Now fit the entropy based weights in your model 
    source(paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/functions/entropy_wt_functions.R"))
    
    # Fitting entropy weights on zip code level variables
    data_ent <- data %>% select(all_of(c("zip", "pm25", confounders)))
    
    # Get rid of repeats since eliminating strata to fit the models
    data_ent <- unique(data_ent)
    
    
    # data_ent_matrix <- model.matrix( ~ -1 + medhouseholdincome + medianhousevalue + poverty + education +
    #                                    pct_owner_occ + year + region, 
    #                                  data = data_ent)
    
    data_ent_matrix <- model.matrix(reformulate(c("-1", confounders)),
                                    data = data_ent)
    
    # data_ent <-
    #   data[, c("medhouseholdincome", "medianhousevalue", "poverty", "education",
    #            "pct_owner_occ", "year", "region", "sex", "race", "dual", "entry_age_break", "followup_year")]
    
    # Create model matrx
    # data_ent_matrix <- model.matrix( ~ -1 + medhouseholdincome + medianhousevalue + poverty + education +
    #                                   pct_owner_occ + year + region + sex + race + dual + entry_age_break + followup_year, 
    #                                 data = data_ent)
    
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
    e <- ebal(data_ent$pm25, c_mat)
    ent_weights <- e$weights
    ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
    
    # # Run through one more iteration to try to get rid of extreme weights
    e <- ebal(data_ent$pm25, c_mat, base_weights = ent_weights)
    ent_weights <- e$weights
    #ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
    
    
    #e <- ebal(data_ent$pm25, c_mat, base_weights = ent_weights)
    #ent_weights <- e$weights
    #ent_weights[ent_weights > quantile(ent_weights, 0.995)] <- quantile(ent_weights, 0.995)
    
    # Now save entropy weights and covariate weights
    cov.wt(cbind(data_ent$pm25, c_mat), ent_weights, TRUE)$cor[, 1]
    cov.wt(cbind(data_ent$pm25, data_ent_matrix), ent_weights, TRUE)$cor[, 1]
    # ent_dataset <- data.table(cbind(pm25_ensemble = data_ent$pm25_ensemble, data_ent_matrix, ent_weights))
    
    # Add the entropy weights to this dataset
    data_ent$ent_weight <- ent_weights
    
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
  
  } 
  
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
                          nthread = 10, 
                          trim = T,
                          data_ent = data_ent,
                          data_confounders = data_confounders)
  
  pseudo_pop <- grid_min_results$pseudo_pop
  
  saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
  
  # pseudo_pop <- 
  #   generate_pseudo_pop(data$log_mort,
  #                       data$pm25_ensemble,
  #                       data.frame(data_ent),
  #                       ci_appr = "matching",
  #                       pred_model = "sl",
  #                       gps_model = "parametric",
  #                       use_cov_transform = TRUE,
  #                       transformers = list("pow2", "pow3", "abs", "scale"),
  #                       trim_quantiles = c(0.025, 0.975),
  #                       optimized_compile = TRUE,
  #                       sl_lib = c("m_xgboost"),
  #                       params = list(xgb_nrounds = c(10, 20, 50, 200),
  #                                     xgb_eta = c(0.025, 0.05, 0.1, 0.2, 0.3),
  #                                     max_depth = c(4, 5, 6)),
  #                       covar_bl_method = "absolute",
  #                       covar_bl_trs = 0.1,
  #                       covar_bl_trs_type = "mean",
  #                       max_attempt = 20,
  #                       matching_fun = "matching_l1",
  #                       delta_n = 0.16,
  #                       scale = 1,
  #                       nthread = 40)
  # 
  # 
  # saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))
  
  #pseudo_pop <- readRDS(file = paste0(out_dir, "pseudo_pop.RDS"))
  
  # Now fit semi-parametric here 
  #out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
  
  # Now add back on zip code?
  psuedo_pop_frame <- data.table(pseudo_pop)

  causalgps_join <- 
    data_ent %>% 
    select(zip, year) %>%
    mutate(causal_weight = psuedo_pop_frame$counter_weight)
  
  data <- 
    data %>% 
    left_join(causalgps_join, by = c("zip", "year")) %>%
    mutate(causal_pop_weight = time_count * causal_weight) # Multiply survey weights
  
  # causal_gps <- 
  #   mgcv::bam(formula = Y ~ s(w, bs = 'cr', k = 4),
  #             family = "gaussian",
  #             data = data,
  #             weights = causal_pop_weight)
  
  
  # causal_gps <-
  #   mgcv::bam(log_mort ~ s(pm25_ensemble, bs = 'cr', k = 4) + medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +
  #               as.factor(year) + as.factor(region) +
  #               as.factor(sex) + as.factor(race) + as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year),
  #           data = data, weights = causal_pop_weight)
  
  # Now fit causal model with strata specific variables in the outcome model
  causal_gps <- 
    mgcv::bam(reformulate(c("s(pm25, bs = 'cr', k = 4)", strata_var, confounders), response = "log_mort"),
              weight = causal_pop_weight, data = data)
  
  
  # causal_gps2 <- 
  #   mgcv::bam(formula = Y ~ s(w, bs = 'ps'),
  #             family = "gaussian",
  #             data = psuedo_pop_frame,
  #             weights = counter_weight)
  # 
  # causal_gps3 <- 
  #   mgcv::bam(formula = Y ~ s(w, bs = 'cr'),
  #             family = "gaussian",
  #             data = psuedo_pop_frame,
  #             weights = counter_weight)
  # 
  # causal_gps4 <- 
  #   mgcv::bam(formula = Y ~ s(w, bs = 'cr', k = 6),
  #             family = "gaussian",
  #             data = psuedo_pop_frame,
  #             weights = counter_weight)
  # 
  # viz_fit1 <- mgcViz::getViz(causal_gps)
  # viz_fit2 <- mgcViz::getViz(causal_gps4)
  # 
  # # Make plotGAM objects
  # trt_fit1 <- plot(viz_fit1, allTerms = T) + l_fitLine()
  # trt_fit2 <- plot(viz_fit2, allTerms = T) + l_fitLine()
  # 
  # exam_dat <- 
  #   bind_rows(trt_fit1[["plots"]][[1]]$data$fit %>% mutate(fit = "Fit 1"), 
  #             trt_fit2[["plots"]][[1]]$data$fit %>% mutate(fit = "Fit 2"))
  # 
  # 
  # ggplot(data = exam_dat, aes(x = x, y = y, colour = fit)) +
  #   geom_line() +
  #   labs(x = "Examination", y = "s(Examination)") + 
  #   theme_bw() +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  saveRDS(causal_gps, file = paste0(out_dir, "causal_fit.RDS"))
  
  message("Done with fitting models")
}

# load the models if doing post processing only
linear_fit <- readRDS(file = paste0(out_dir, "linear.RDS"))
gam_fit <- readRDS(file = paste0(out_dir, "gam.RDS"))
linear_ent <- readRDS(file = paste0(out_dir, "linear_ent.RDS"))
gam_ent <- readRDS(file = paste0(out_dir, "gam_ent.RDS"))

# Load causal fit
pseudo_pop <- readRDS(file = paste0(out_dir, "pseudo_pop.RDS"))
causal_gps <- readRDS(file = paste0(out_dir, "causal_fit.RDS"))

# only plot from first to 99th percentile in accordance with past papers
# This corresponds to 

# min_pm <- quantile(data_ent$pm25, 0.01)
# max_pm <- quantile(data_ent$pm25, 0.99)

# plot entire curve
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

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

data_prediction <- readRDS(file = paste0(out_dir, "data_prediction.RDS"))

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

# Plot the contrasts
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

# Add some processing to the data, do I want to run on age-adjusted mortality rate? Or log of mortality rate?

# I think I need zip code population per year (at the zip code level?)

# fit_causal_gps <- function(data) {
#   
#   # Now fit most optimized causal pathway
#   # Set grid to tune over 
#   nrounds = 300
#   tune_grid <- expand.grid(
#     nrounds = nrounds,
#     eta = c(0.025, 0.05, 0.1),
#     max_depth = c(2, 3),
#     delta = seq(0.2, 2, by = 0.2)
#   )
#   
#   # Make list to pass to parLapply
#   tune_grid_list <- as.list(as.data.frame(t(tune_grid)))
#   
#   # Wrapper function for running causalGPS
#   wrapper_func <- function(tune_grid){
#     pseudo_pop <- generate_pseudo_pop(data$log_mort,
#                                       data$pm25_ensemble,
#                                       data.frame(dplyr::select(data, medhouseholdincome, medianhousevalue, 
#                                                     poverty, education, pct_owner_occ,
#                                                     year,region, sex, race, dual, entry_age_break)),
#                                       ci_appr = "matching",
#                                       pred_model = "sl",
#                                       gps_model = "parametric",
#                                       use_cov_transform = FALSE,
#                                       transformers = list("pow2", "pow3", "abs", "scale"),
#                                       trim_quantiles = c(0.025,0.975),
#                                       optimized_compile = T,
#                                       sl_lib = c("m_xgboost"),
#                                       params = list(xgb_rounds = tune_grid[[1]],
#                                                     xgb_eta = tune_grid[[2]],
#                                                     xgb_max_depth = tune_grid[[3]]),
#                                       covar_bl_method = "absolute",
#                                       covar_bl_trs = 0.1,
#                                       covar_bl_trs_type = "mean",
#                                       max_attempt = 1,
#                                       matching_fun = "matching_l1",
#                                       delta_n = tune_grid[[4]],
#                                       scale = 1,
#                                       nthread = 2)
#     
#     
#     results <- 
#       data.frame(nrounds = tune_grid[[1]], 
#                  eta = tune_grid[[2]],
#                  max_depth = tune_grid[[3]],
#                  delta = tune_grid[[4]],
#                  scale = 1,
#                  mean_corr = pseudo_pop$adjusted_corr_results$mean_absolute_corr,
#                  max_corr = pseudo_pop$adjusted_corr_results$maximal_absolute_corr)
#     return(list(results, pseudo_pop))
#     
#   }
#   
#   # Run these wrapper functions in parallel
#   # Consider running optimized vs non optimized version through simulations? Could try to get that working today 
#   
#   # currently not running parallel so commenting this out 
#   # cl <- parallel::makeCluster(10, type="PSOCK")
#   # parallel::clusterExport(cl=cl,
#   #                         varlist = c("generate_pseudo_pop",
#   #                                     "wrapper_func",
#   #                                     "data_ipw"
#   #                         ),
#   #                         envir=environment())
#   # 
#   # pseudo_pop_list  <- parallel::parLapply(cl,tune_grid_list, wrapper_func)
#   # parallel::stopCluster(cl)
#   
#   pseudo_pop_list <- mclapply(tune_grid_list, wrapper_func, mc.cores = 4)
#   
#   # If you don't want to run in parallel, use this 
#   #pseudo_pop_list <- mclapply(tune_grid_list, mc.cores = 4, wrapper_func)
#   corr_search <- do.call("rbind", (lapply(pseudo_pop_list, function(x) x[[1]])))
#   
#   # Extract minimum as your result
#   pseudo_pop_tuned <- pseudo_pop_list[[which.min(corr_search$mean_corr)]][[2]]
#   
#   # Now fit semi-parametric here 
#   causal_gps_tuned <- gam::gam(formula = Y ~ s(w, df = 3), family = "gaussian", data = data.frame(pseudo_pop_tuned$pseudo_pop), weights = counter)
#   
# }
