rm(list = ls())

# Set library path


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
only_causalgps = as.logical(args[[3]])
delta_n <- as.numeric(args[[4]])

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/input_data.RDS"))

data[, year := factor(year)]
data[, followup_year := factor(followup_year)]

# Create output directory for models
out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
dir.create(out_dir)

# Fit a glm
if (!only_causalgps) {
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
  
  saveRDS(gam_fit, file = paste0(out_dir, "gam.RDS"))
  
  # gam_fit_bs <- 
  #   mgcv::bam(log_mort ~ s(pm25_ensemble, bs = 'bs', k = 20) + medhouseholdincome + medianhousevalue + poverty + education + pct_owner_occ +as.factor(year) + as.factor(region) + 
  #        as.factor(sex) + as.factor(race) + as.factor(dual) + as.factor(entry_age_break) + as.factor(followup_year),
  #      data = data)
  # 
  # saveRDS(gam_fit_bs, file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/gam_bs.RDS"))
  # 
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
  
  data_ent <-
    data[, c("medhouseholdincome", "medianhousevalue", "poverty", "education",
             "pct_owner_occ", "year", "region", "sex", "race", "dual", "entry_age_break", "followup_year")]
  
  # Create model matrx
  data_ent_matrix <- model.matrix( ~ -1 + medhouseholdincome + medianhousevalue + poverty + education +
                                    pct_owner_occ + year + region + sex + race + dual + entry_age_break + followup_year, 
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
  saveRDS(ent_weights, file = paste0(out_dir, "entropy_weights_raw.RDS"))
   
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

} 

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
                      covar_bl_trs = 0.5,
                      covar_bl_trs_type = "mean",
                      max_attempt = 1,
                      matching_fun = "matching_l1",
                      delta_n = 0.16,
                      scale = 1,
                      nthread = 40)

saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop.RDS"))

# Now create year in entropy weights
data_ent[, year := as.numeric(as.character((year)))]

pseudo_pop <- 
  generate_pseudo_pop(data$log_mort,
                      data$pm25_ensemble,
                      data.frame(data_ent),
                      ci_appr = "matching",
                      pred_model = "sl",
                      gps_model = "parametric",
                      use_cov_transform = TRUE,
                      transformers = list("pow2", "pow3", "abs", "scale"),
                      trim_quantiles = c(0.01, 0.99),
                      optimized_compile = TRUE,
                      sl_lib = c("m_xgboost"),
                      params = list(xgb_nrounds = c(50),
                                    xgb_eta = c(.3),
                                    max_depth = c(6),
                                    "xgb_min_child_weight" = 1),
                      covar_bl_method = "absolute",
                      covar_bl_trs = 0.3,
                      covar_bl_trs_type = "mean",
                      max_attempt = 1,
                      matching_fun = "matching_l1",
                      delta_n = 0.16,
                      scale = 1,
                      nthread = 40)

saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop2.RDS"))

# Now run again with defaults
pseudo_pop <- 
  generate_pseudo_pop(data$log_mort,
                      data$pm25_ensemble,
                      data.frame(data_ent),
                      ci_appr = "matching",
                      pred_model = "sl",
                      gps_model = "parametric",
                      use_cov_transform = TRUE,
                      transformers = list("pow2", "pow3", "abs", "scale"),
                      trim_quantiles = c(0.025, 0.975),
                      optimized_compile = TRUE,
                      sl_lib = c("m_xgboost"),
                      params = list(xgb_nrounds = c(10, 20, 50, 200),
                                    xgb_eta = c(0.025, 0.05, 0.1, 0.2, 0.3),
                                    max_depth = c(4, 5, 6)),
                      covar_bl_method = "absolute",
                      covar_bl_trs = 0.1,
                      covar_bl_trs_type = "mean",
                      max_attempt = 20,
                      matching_fun = "matching_l1",
                      delta_n = 0.16,
                      scale = 1,
                      nthread = 40)


saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop3.RDS"))

# Now run with trimmed
pseudo_pop <- 
  generate_pseudo_pop(data$log_mort,
                      data$pm25_ensemble,
                      data.frame(data_ent),
                      ci_appr = "matching",
                      pred_model = "sl",
                      gps_model = "parametric",
                      use_cov_transform = TRUE,
                      transformers = list("pow2", "pow3", "abs", "scale"),
                      trim_quantiles = c(0.01, 0.99),
                      optimized_compile = TRUE,
                      sl_lib = c("m_xgboost"),
                      params = list(xgb_nrounds = c(10, 20, 50, 200),
                                    xgb_eta = c(0.025, 0.05, 0.1, 0.2, 0.3),
                                    max_depth = c(4, 5, 6)),
                      covar_bl_method = "absolute",
                      covar_bl_trs = 0.1,
                      covar_bl_trs_type = "mean",
                      max_attempt = 20,
                      matching_fun = "matching_l1",
                      delta_n = 0.16,
                      scale = 1,
                      nthread = 40)


saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop4.RDS"))

# Now turn into factor again
data_ent[, year := as.factor(year)]

pseudo_pop <- 
  generate_pseudo_pop(data$log_mort,
                      data$pm25_ensemble,
                      data.frame(data_ent),
                      ci_appr = "matching",
                      pred_model = "sl",
                      gps_model = "parametric",
                      use_cov_transform = TRUE,
                      transformers = list("pow2", "pow3", "abs", "scale"),
                      trim_quantiles = c(0.025, 0.975),
                      optimized_compile = TRUE,
                      sl_lib = c("m_xgboost"),
                      params = list(xgb_nrounds = c(10, 20, 50, 200),
                                    xgb_eta = c(0.025, 0.05, 0.1, 0.2, 0.3),
                                    max_depth = c(4, 5, 6)),
                      covar_bl_method = "absolute",
                      covar_bl_trs = 0.1,
                      covar_bl_trs_type = "mean",
                      max_attempt = 20,
                      matching_fun = "matching_l1",
                      delta_n = 0.16,
                      scale = 1,
                      nthread = 40)


saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop5.RDS"))

# Different trim
pseudo_pop <- 
  generate_pseudo_pop(data$log_mort,
                      data$pm25_ensemble,
                      data.frame(data_ent),
                      ci_appr = "matching",
                      pred_model = "sl",
                      gps_model = "parametric",
                      use_cov_transform = TRUE,
                      transformers = list("pow2", "pow3", "abs", "scale"),
                      trim_quantiles = c(0.01, 0.99),
                      optimized_compile = TRUE,
                      sl_lib = c("m_xgboost"),
                      params = list(xgb_nrounds = c(10, 20, 50, 200),
                                    xgb_eta = c(0.025, 0.05, 0.1, 0.2, 0.3),
                                    max_depth = c(4, 5, 6)),
                      covar_bl_method = "absolute",
                      covar_bl_trs = 0.1,
                      covar_bl_trs_type = "mean",
                      max_attempt = 20,
                      matching_fun = "matching_l1",
                      delta_n = 0.16,
                      scale = 1,
                      nthread = 40)


saveRDS(pseudo_pop, file = paste0(out_dir, "pseudo_pop6.RDS"))

#pseudo_pop <- readRDS(file = paste0(out_dir, "pseudo_pop.RDS"))

# Now fit semi-parametric here 
#out_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
psuedo_pop_frame <- data.table(pseudo_pop$pseudo_pop)

causal_gps <- 
  mgcv::bam(formula = Y ~ s(w, bs = 'cr', k = 4),
            family = "gaussian",
            data = psuedo_pop_frame,
            weights = counter_weight)

saveRDS(causal_gps, file = paste0(out_dir, "causal_fit.RDS"))

message("Done with fitting models")

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
