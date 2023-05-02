# This script is for deploying a grid search to find the right hyperparameter for xgboost settings for my simulation 
rm(list = ls())
library(data.table)
library(tidyverse)
library(caret)
library(CausalGPS)

source("~/Desktop/Francesca_research/Simulation_studies/simulation_functions.R")
cov_function <- function(confounders) as.vector(-0.8 + matrix(c(0.1, 0.1, -0.1, 0.2, 0.1, 0.1), nrow = 1) %*% t(confounders))
scale_exposure <- function(x){20 * (x-min(x))/(max(x)-min(x))}

sample_size = 1000
cf <- mvrnorm(n = sample_size,
              mu = rep(0, 4),
              Sigma = diag(4))
cf5 <- sample(c((-2):2), sample_size, replace = T)
cf6 <- runif(sample_size, min = -3, max = 3)
confounders = cbind(cf, cf5, cf6)
colnames(confounders) = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")

# Now doing it for first exposure
exposure = scale_exposure(cov_function(confounders)) + rnorm(sample_size, mean = 0, sd = sqrt(10))

# Now generate a dataset (must set these to be correct)
data_example <- data_generate_a(sample_size = sample_size, exposure = exposure, confounders = confounders, 
                                exposure_relationship = "linear", outcome_relationship = "linear",
                                family = "gaussian", confounder_mult = 1)
data_gps = data_example %>% dplyr::select(-Y)

# Method to tune XGBoost
tune_control <-
  caret::trainControl(method = "cv",
                      number = 5,
                      verboseIter = FALSE,
                      allowParallel = TRUE)

nrounds = 700
tune_grid <- expand.grid(
  nrounds = seq(from = 200, to = nrounds, by = 50),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(2, 3, 4),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

# Initial test of tuning parameters
xgb_tune <-
  train(exposure ~ ., data = data_gps, method = "xgbTree",
        trControl = tune_control, tuneGrid = tune_grid, verbosity = 0)

# Tuning parameters we will use to search for optimal delta 
xgb_tune$bestTune

corr_search_results <- list()
for (j in 1:10) {
  cf <- mvrnorm(n = sample_size,
                mu = rep(0, 4),
                Sigma = diag(4))
  cf5 <- sample(c((-2):2), sample_size, replace = T)
  cf6 <- runif(sample_size, min = -3, max = 3)
  confounders = cbind(cf, cf5, cf6)
  colnames(confounders) = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")
  
  # Now doing it for first exposure
  scale_exposure <- function(x){20 * (x-min(x))/(max(x)-min(x))}
  exposure_1 = scale_exposure(cov_function(confounders)) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  exposure_2 = scale_exposure(cov_function(confounders)) + (sqrt(5)) * rt(sample_size, df = 3)
  exposure_3 = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  exposure_4 = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2 + 0.5 * (confounders[, "cf1"]) * (confounders[, "cf5"])) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  
  corr_search_delta <- 
    rbindlist(lapply(c(1:4), function(exp) {
      data_example <- data_generate_a(sample_size = sample_size, exposure = get(paste0("exposure_", exp)), confounders = confounders, 
                                      exposure_relationship = "linear", outcome_relationship = "linear",
                                      family = "gaussian", confounder_mult = 1)
      
      delta_search <- seq(1, 4, by = 0.2)
      # Maybe start with delta optimizing first, then do the xgboost grid search
      corr_search_delta <- 
        rbindlist(mclapply(1:length(delta_search), mc.cores = 2, function(i){
          pseudo_pop <- generate_pseudo_pop(data_example$Y,
                                            data_example$exposure,
                                            data.frame(data_example[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                            ci_appr = "matching",
                                            pred_model = "sl",
                                            gps_model = "non-parametric",
                                            use_cov_transform = FALSE,
                                            transformers = list("pow2", "pow3", "abs", "scale"),
                                            trim_quantiles = c(0.025,0.975),
                                            optimized_compile = F,
                                            sl_lib = c("m_xgboost"),
                                            params = list(xgb_max_depth = xgb_tune$bestTune$max_depth,
                                                          xgb_rounds = xgb_tune$bestTune$nrounds,
                                                          xgb_eta = xgb_tune$bestTune$eta),
                                            covar_bl_method = "absolute",
                                            covar_bl_trs = 0.1,
                                            covar_bl_trs_type = "mean",
                                            max_attempt = 1,
                                            matching_fun = "matching_l1",
                                            delta_n = delta_search[i],
                                            scale = 1,
                                            nthread = 1)
          
          results <- 
            data.table(delta = delta_search[i],
                       mean_corr = pseudo_pop$adjusted_corr_results$mean_absolute_corr,
                       max_corr = pseudo_pop$adjusted_corr_results$maximal_absolute_corr)
          return(results)
        }))
      corr_search_delta$exposure_type = exp
      return(corr_search_delta)
  }))
  corr_search_delta$sim <- j
  corr_search_results[[j]] <- corr_search_delta
}

# Save the 
corr_search_results %>% 
  rbindlist() %>% 
  ggplot() + 
  geom_point(aes(x = delta, y = mean_corr)) + 
  geom_smooth(aes(x = delta, y = mean_corr)) + 
  facet_wrap(~ exposure_type)
  
# Plot of correlation by delta 
corr_search_delta %>% 
  ggplot() + 
  geom_point(aes(x = delta, y = mean_corr)) + 
  geom_smooth(aes(x = delta, y = mean_corr))


min_delta <- corr_search_delta %>% filter(mean_corr == min(mean_corr)) %>% pull(delta)


# Now varying the xgboost parameters
corr_search_results_xgb <- list()
for (j in 3:10) {
  cf <- mvrnorm(n = sample_size,
                mu = rep(0, 4),
                Sigma = diag(4))
  cf5 <- sample(c((-2):2), sample_size, replace = T)
  cf6 <- runif(sample_size, min = -3, max = 3)
  confounders = cbind(cf, cf5, cf6)
  colnames(confounders) = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")
  
  # Now doing it for first exposure
  scale_exposure <- function(x){20 * (x-min(x))/(max(x)-min(x))}
  exposure_1 = scale_exposure(cov_function(confounders)) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  exposure_2 = scale_exposure(cov_function(confounders)) + (sqrt(5)) * rt(sample_size, df = 3)
  exposure_3 = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  exposure_4 = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2 + 0.5 * (confounders[, "cf1"]) * (confounders[, "cf5"])) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  
  corr_search_xgb_outer <- 
    rbindlist(lapply(c(1:4), function(exp) {
      data_example <- data_generate_a(sample_size = sample_size, exposure = get(paste0("exposure_", exp)), confounders = confounders, 
                                      exposure_relationship = "linear", outcome_relationship = "linear",
                                      family = "gaussian", confounder_mult = 1)
      
      nrounds = 300
      tune_grid <- expand.grid(
        nrounds = nrounds,
        eta = c(0.025, 0.05, 0.1),
        max_depth = c(2, 3),
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
      )
      
      # noticed the number of rounds is less important than the eta and learning rate
      
      #delta_search <- seq(1, 4, by = 0.1)
      # Maybe start with delta optimizing first, then do the xgboost grid search
      corr_search_xgb <-
        rbindlist(apply(tune_grid, 1, function(i){
          pseudo_pop <- generate_pseudo_pop(data_example$Y,
                                            data_example$exposure,
                                            data.frame(data_example[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                            ci_appr = "matching",
                                            pred_model = "sl",
                                            gps_model = "non-parametric",
                                            use_cov_transform = FALSE,
                                            transformers = list("pow2", "pow3", "abs", "scale"),
                                            trim_quantiles = c(0.025,0.975),
                                            optimized_compile = F,
                                            sl_lib = c("m_xgboost"),
                                            params = list(xgb_max_depth = i[[3]],
                                                          xgb_rounds = i[[1]],
                                                          xgb_eta = i[[2]]),
                                            covar_bl_method = "absolute",
                                            covar_bl_trs = 0.1,
                                            covar_bl_trs_type = "mean",
                                            max_attempt = 1,
                                            matching_fun = "matching_l1",
                                            delta_n = 1,
                                            scale = 1,
                                            nthread = 1)
          
          results <- 
            data.table(max_depth = i[[3]],
                       nrounds = i[[1]], 
                       eta = i[[2]], 
                       mean_corr = pseudo_pop$adjusted_corr_results$mean_absolute_corr,
                       max_corr = pseudo_pop$adjusted_corr_results$maximal_absolute_corr)
          return(results)
        }))
      corr_search_xgb$exposure_type = exp
      return(corr_search_xgb)
    }))
  corr_search_xgb_outer$sim <- j
  corr_search_results_xgb[[j]] <- corr_search_xgb_outer
}


# Save the 
corr_search_results_xgb %>% 
  rbindlist() %>% 
  ggplot() + 
  geom_point(aes(x = max_depth, y = mean_corr, color = factor(max_depth))) + 
  geom_smooth(aes(x = max_depth, y = mean_corr), method = "lm") + 
  facet_wrap(~ exposure_type)


# Make sure it can reach appropriate 0.1 ----------
# Just an experiment
sample_size = 1000
#corr_search_results_all <- list()
for (j in 4:10) {
  cf <- mvrnorm(n = sample_size,
                mu = rep(0, 4),
                Sigma = diag(4))
  cf5 <- sample(c((-2):2), sample_size, replace = T)
  cf6 <- runif(sample_size, min = -3, max = 3)
  confounders = cbind(cf, cf5, cf6)
  colnames(confounders) = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")
  
  # Now doing it for first exposure
  scale_exposure <- function(x){20 * (x-min(x))/(max(x)-min(x))}
  exposure_1 = scale_exposure(cov_function(confounders)) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  exposure_2 = scale_exposure(cov_function(confounders)) + (sqrt(5)) * rt(sample_size, df = 3)
  exposure_3 = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  exposure_4 = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2 + 0.5 * (confounders[, "cf1"]) * (confounders[, "cf5"])) + rnorm(sample_size, mean = 0, sd = sqrt(10))
  
  corr_search_all_outer <- 
    rbindlist(lapply(c(1:4), function(exp) {
      data_example <- data_generate_a(sample_size = sample_size, exposure = get(paste0("exposure_", exp)), confounders = confounders, 
                                      exposure_relationship = "linear", outcome_relationship = "linear",
                                      family = "gaussian", confounder_mult = 1)
      
      nrounds = 300
      tune_grid <- expand.grid(
        nrounds = nrounds,
        eta = c(0.025, 0.05, 0.1),
        max_depth = c(2, 3),
        delta = seq(0.2, 2, by = 0.2)
      )
      
      tune_grid_list <- as.list(as.data.frame(t(tune_grid)))
      
      # noticed the number of rounds is less important than the eta and learning rate, therefore keep that fixed at 300
      wrapper_func <- function(tune_grid){
        pseudo_pop <- generate_pseudo_pop(data_example$Y,
                                          data_example$exposure,
                                          data.frame(data_example[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                          ci_appr = "matching",
                                          pred_model = "sl",
                                          gps_model = "non-parametric",
                                          use_cov_transform = FALSE,
                                          transformers = list("pow2", "pow3", "abs", "scale"),
                                          trim_quantiles = c(0.025,0.975),
                                          optimized_compile = T,
                                          sl_lib = c("m_xgboost"),
                                          params = list(xgb_rounds = tune_grid[[1]],
                                                        xgb_eta = tune_grid[[2]],
                                                        xgb_max_depth = tune_grid[[3]]),
                                          covar_bl_method = "absolute",
                                          covar_bl_trs = 0.1,
                                          covar_bl_trs_type = "mean",
                                          max_attempt = 1,
                                          matching_fun = "matching_l1",
                                          delta_n = tune_grid[[4]],
                                          scale = 1,
                                          nthread = 1)
        
        
        results <- 
          data.frame(nrounds = tune_grid[[1]], 
                     eta = tune_grid[[2]],
                     max_depth = tune_grid[[3]],
                     delat = tune_grid[[4]],
                     lambda = 1,
                     mean_corr = pseudo_pop$adjusted_corr_results$mean_absolute_corr,
                     max_corr = pseudo_pop$adjusted_corr_results$maximal_absolute_corr)
        return(results)
        
      }
      
      # Run these wrapper functions in parallel
      # Consider running optimized vs non optimized version through simulations? Could try to get that working today 
      cl <- parallel::makeCluster(12, type="PSOCK")
      parallel::clusterExport(cl=cl,
                              varlist = c("generate_pseudo_pop",
                                          "wrapper_func",
                                          "data_example"
                              ),
                              envir=environment())
      
      pseudo_pop_list_2  <- parallel::parLapply(cl,tune_grid_list, wrapper_func)
      parallel::stopCluster(cl)
      corr_search <- do.call("rbind", pseudo_pop_list_2)
      corr_search$exposure <- exp
      return(corr_search)
    }))
  corr_search_all_outer$sim <- j
  corr_search_results_all[[j]] <- corr_search_all_outer
}


# Now implement this in the GPS package!!! That is the next step 
      

# now extract what you want to extract!! Also check the memory now!
      
      corr_search_all <-
        rbindlist(apply(tune_grid, 1, function(i){
          pseudo_pop <- generate_pseudo_pop(data_example$Y,
                                            data_example$exposure,
                                            data.frame(data_example[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                            ci_appr = "matching",
                                            pred_model = "sl",
                                            gps_model = "non-parametric",
                                            use_cov_transform = FALSE,
                                            transformers = list("pow2", "pow3", "abs", "scale"),
                                            trim_quantiles = c(0.025,0.975),
                                            optimized_compile = T,
                                            sl_lib = c("m_xgboost"),
                                            params = list(xgb_max_depth = i[[3]],
                                                          xgb_rounds = i[[1]],
                                                          xgb_eta = i[[2]]),
                                            covar_bl_method = "absolute",
                                            covar_bl_trs = 0.1,
                                            covar_bl_trs_type = "mean",
                                            max_attempt = 1,
                                            matching_fun = "matching_l1",
                                            delta_n = i[[4]],
                                            scale = 1,
                                            nthread = 1)
          
          results <- 
            data.table(max_depth = i[[3]],
                       nrounds = i[[1]], 
                       eta = i[[2]], 
                       delat = i[[4]],
                       lambda = 1,
                       mean_corr = pseudo_pop$adjusted_corr_results$mean_absolute_corr,
                       max_corr = pseudo_pop$adjusted_corr_results$maximal_absolute_corr)
          return(results)
        }))
      corr_search_all$exposure_type = exp
      return(corr_search_all)
    }))
  corr_search_all_outer$sim <- j
  corr_search_results_all[[j]] <- corr_search_all_outer
}

# Save this version to talk about with Francesca
saveRDS(corr_search_results_all, file = "~/Desktop/Francesca_research/Simulation_studies/corr_results.R")

corr_search_results_all <- readRDS(file = "~/Desktop/Francesca_research/Simulation_studies/corr_results.R")
# Save the 
corr_search_results_all %>%
  rbindlist() %>% 
  ggplot() + 
  geom_point(aes(x = lambda, y = mean_corr)) + 
  geom_smooth(aes(x = lambda, y = mean_corr), method = "lm") + 
  facet_wrap(~ exposure_type)

corr_search_results_all %>% 
  rbindlist() %>% 
  filter(exposure_type == 1) %>% 
  lm(mean_corr ~ eta + delat + lambda + max_depth, data = .)


# In case I want to do a final tune (right now not interested)
corr_search_delta %>% filter(mean_corr < quantile(corr_search_delta$mean_corr, 0.25))




# Ok try to runn this in parallel

delta_n <- c(1:10)*0.1
m_d <- generate_syn_data(sample_size = 200)

wrapper_func <- function(delta_n){
  pseudo_pop <- generate_pseudo_pop(m_d$Y,
                                    m_d$treat,
                                    m_d[c("cf1","cf2","cf3","cf4","cf5","cf6")],
                                    ci_appr = "matching",
                                    pred_model = "sl",
                                    sl_lib = c("m_xgboost"),
                                    params = list(xgb_nrounds=c(10,20,30),
                                                  xgb_eta=c(0.1,0.2,0.3)),
                                    nthread = 1,
                                    covar_bl_method = "absolute",
                                    covar_bl_trs = 0.1,
                                    covar_bl_trs_type = "mean",
                                    max_attempt = 1,
                                    matching_fun = "matching_l1",
                                    delta_n = delta_n,
                                    scale = 0.5)
  
  return(pseudo_pop)
}

# Function to collect the results
print_mean_covariate_balance <- function(delta_n, object_list){
  m_cov <- unlist(lapply(object_list,
                         function(x){x$adjusted_corr_results$mean_absolute_corr}))
  results <- data.frame(delta_n, m_cov)
  print(results)
}


# with lapply ------------------------------------------------------------------
pseudo_pop_list_1  <- lapply(delta_n, wrapper_func)

# mean covariate balance
print_mean_covariate_balance(delta_n = delta_n, object_list = pseudo_pop_list_1)


# with parLapply ---------------------------------------------------------------
cl <- parallel::makeCluster(10, type="PSOCK")
parallel::clusterExport(cl=cl,
                        varlist = c("generate_pseudo_pop",
                                    "wrapper_func",
                                    "m_d"
                        ),
                        envir=environment())

pseudo_pop_list_2  <- parallel::parLapply(cl,delta_n, wrapper_func)
parallel::stopCluster(cl)

# mean covariate balance
print_mean_covariate_balance(delta_n = delta_n, object_list = pseudo_pop_list_2)






# Now select some hyperparameters to iterate over (those that seem to predict the propensity score reasonably well)
# Widdle this down to perform grid search after seeing what is happening
grid_trial_xgboost <- xgb_tune$results %>% filter(RMSE < quantile(xgb_tune$results[["RMSE"]], 0.25)) %>% dplyr::select(eta, max_depth, nrounds) %>% unique()
unique_values_xgboost <- grid_trial_xgboost %>% dplyr::select(eta, max_depth) %>% unique()
grid_trial_xgboost_real <- 
  rbindlist(apply(unique_values_xgboost, 1, function(x){
  grid_trial_xgboost %>% filter(eta == x[[1]], max_depth == x[[2]]) %>% .[1,]
}))

grid_trial <- expand.grid.df(data.frame(grid_trial_xgboost_real), data.frame(delta = seq(0.9, 2, by = 0.1)))

corr_grid_search <- 
  rbindlist(lapply(1:nrow(grid_trial), function(i){
    pop <- generate_pseudo_pop(data_example$Y,
                                      data_example$exposure,
                                      data.frame(data_example[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                      ci_appr = "matching",
                                      pred_model = "sl",
                                      gps_model = "non-parametric",
                                      use_cov_transform = FALSE,
                                      transformers = list("pow2", "pow3", "abs", "scale"),
                                      trim_quantiles = c(0.05,0.95),
                                      optimized_compile = F,
                                      sl_lib = c("m_xgboost"),
                                      params = list(xgb_max_depth = grid_trial$max_depth[i],
                                                    xgb_rounds = grid_trial$nrounds[i],
                                                    xgb_eta = grid_trial$eta[i]),
                                      covar_bl_method = "absolute",
                                      covar_bl_trs = 0.1,
                                      covar_bl_trs_type = "mean",
                                      max_attempt = 1,
                                      matching_fun = "matching_l1",
                                      delta_n = min_delta,
                                      scale = 1,
                                      nthread = 1)
    
    results <- 
      data.table(xgb_max_depth = grid_trial$max_depth[i],
                 xgb_rounds = grid_trial$nrounds[i],
                 xgb_eta = grid_trial$eta[i],
                 mean_corr = pop$adjusted_corr_results$mean_absolute_corr,
                 max_corr = pop$adjusted_corr_results$maximal_absolute_corr)
    return(results)
  }))

corr_grid_search %>% filter(mean_corr == min(mean_corr)) %>% filter(max_corr == min(max_corr))



cutoff <- corr_grid_search$mean_corr %>% quantile(0.2)
final_xgboost_params <- corr_grid_search %>% filter(mean_corr < cutoff) %>% dplyr::select(xgb_max_depth, xgb_rounds, xgb_eta) %>% data.frame()
final_delta_search <- data.frame(delta = seq(0.8, 2.4, by = 0.2))

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
grid_trial2 <- expand.grid.df(final_xgboost_params, final_delta_search)


tune_grid2 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 50),
  eta = xgb_tune$bestTune$eta,
  max_depth = ifelse(xgb_tune$bestTune$max_depth == 2,
                     c(xgb_tune$bestTune$max_depth:4),
                     xgb_tune$bestTune$max_depth - 1:xgb_tune$bestTune$max_depth + 1),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(1, 2, 3),
  subsample = 1
)
xgb_tune2 <-
  train(exposure ~ ., data = data_gps, method = "xgbTree",
        trControl = tune_control, tuneGrid = tune_grid2, verbosity = 0)



tune_grid3 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 50),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = 0,
  colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = c(0.5, 0.75, 1.0)
)

xgb_tune3 <-
  train(exposure ~ ., data = data_gps, method = "xgbTree",
        trControl = tune_control, tuneGrid = tune_grid3, verbosity = 0)




tune_grid4 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 50),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0),
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune4 <-
  train(exposure ~ ., data = data_gps, method = "xgbTree",
        trControl = tune_control, tuneGrid = tune_grid4, verbosity = 0)

tune_grid5 <- expand.grid(
  nrounds = seq(from = 100, to = 10000, by = 100),
  eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = xgb_tune4$bestTune$gamma,
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune5 <-
  train(exposure ~ ., data = data_gps, method = "xgbTree",
        trControl = tune_control, tuneGrid = tune_grid5, verbosity = 0)

plot(xgb_tune5)




xgbt_model$bestTune

predVals <- extractPrediction(list(xgb_tune3))
plotObsVsPred(predVals)





med_pred_gps <- 
  lapply(c(1:4), function(gps_mod) {
    if (gps_mod == 1) {
      exposure = scale_exposure(cov_function(confounders)) + rnorm(sample_size, mean = 0, sd = sqrt(10))
    } else if (gps_mod == 2) {
      exposure = scale_exposure(cov_function(confounders)) + (sqrt(5)) * rt(sample_size, df = 3)
    } else if (gps_mod == 3) {
      exposure = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2) + rnorm(sample_size, mean = 0, sd = sqrt(10))
    } else if (gps_mod == 4) {
      exposure = scale_exposure(cov_function(confounders) + 0.5 * (confounders[, "cf3"]) ^ 2 + 0.5 * (confounders[, "cf1"]) * (confounders[, "cf5"])) + rnorm(sample_size, mean = 0, sd = sqrt(10))
    }