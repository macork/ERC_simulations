# Code used for simulation functions

# Function to generate data
data_generate_a <-
  function(sample_size = 1000, 
           exposure = NA, 
           confounders = NA, 
           exposure_relationship = NA, 
           outcome_relationship = NA, 
           family = NA, 
           outcome_sd = 10) {
    
    # Abbreviate for clarity 
    cf <- confounders
    
    # For now mostly using gaussian family, though poisson has been used in the past as well
    if (family == "gaussian") {
      
      # Add on exposure effect depending on relationship
      if (exposure_relationship == "linear") {
        transformed_exp <- exposure
        Y = transformed_exp
      } else if (exposure_relationship == "sublinear") {
        transformed_exp = 5 * log(exposure + 1)
        Y = transformed_exp
      } else if (exposure_relationship == "threshold") {
        # Set threshold at 5
        transformed_exp <- exposure 
        transformed_exp[transformed_exp <= 5] <- 0
        transformed_exp[transformed_exp > 5] <- transformed_exp[transformed_exp > 5] - 5
        transformed_exp <- 1.5 * transformed_exp
        Y = transformed_exp
      }
      
      # Add on effect of the confounders
      if (outcome_relationship == "linear") {
        Y = Y + as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(cf) + rnorm(sample_size, mean = 0, sd = outcome_sd))
      } else if (outcome_relationship == "interaction") {
        # problem with the threshold exposure scenario happening 
        Y = Y + as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(cf) - 
                             transformed_exp*(-0.1 * cf[, 1] + 0.1 * cf[, 3]^2 + 0.1*cf[, 4] + 0.1*cf[, 5]) + rnorm(sample_size, mean = 0, sd = outcome_sd))
      } else {
        stop("Outcome relationship for now is only between linear and interaction")
      }
      
    } else if (family == "poisson") {
      # linear case
      lambdas = exp(2 +  as.vector(matrix(c(0.2, 0.2, 0.3, -0.1, -0.2, 0.2), nrow = 1) %*% t(cf)) + 0.1 * exposure)
      Y = rpois(n = 1000, lambda = lambdas)
      
      # sublinear case 
      lambdas = exp(as.numeric(unlist(2 + 0.2 * cf[, 1] + 0.2 * cf[, 2] + 0.3 * cf[, 3] -0.1 * cf[, 4] - 0.2 * cf[, 5] + 0.2 * cf[, 6] + 1.5 * log10(exposure + 1))))
      Y = rpois(n = 1000, lambda = lambdas)
      
      # threshold case
      thresh_exp <- exposure 
      thresh_exp[thresh_exp <= 5] <- 0
      thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
      lambdas = exp(as.numeric(unlist(2 + 0.2 * cf[, 1] + 0.2 * cf[, 2] + 0.3 * cf[, 3] -0.1 * cf[, 4] - 0.2 * cf[, 5] + 0.2 * cf[, 6] + 0.1 * thresh_exp)))
      Y = rpois(n = 1000, lambda = lambdas)
    }
    
    # Save into one simulated dataset
    simulated_data <- data.table(cbind(Y, exposure, cf))
    colnames(simulated_data)[3:8] <- c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")
    return(simulated_data)
  }

# This function fits all estimators to the data and extracts the absolute bias, RMSE
metrics_from_data <- function(exposure = NA, 
                              exposure_relationship = "linear", # Underlying exposure-response relationship
                              outcome_relationship = "linear", #linear or interaction
                              sample_size = 1000, 
                              confounders = NA,
                              family = "gaussian", # Family of simulation (poisson or gaussian)
                              eschif_draws = NULL, # Should eSCHIF be used in model run
                              adjust_confounder = T, # Should confounders be included in outcome model
                              causal_gps = F) {
  
  # Simulate data given exposure, confounders and sample size
  data_example <- data_generate_a(sample_size = sample_size, 
                                  exposure = exposure, 
                                  confounders = confounders, 
                                  exposure_relationship = exposure_relationship, 
                                  outcome_relationship = outcome_relationship,
                                  family = family)
  
  # Now fit a linear model
  if (adjust_confounder == T) {
    linear_fit <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                      data = data_example, 
                      family = family)
  } else {
    linear_fit <- glm(Y ~ exposure, data = data_example, family = family)
  }
  
  # Now fit a GAM model 
  if (adjust_confounder == T) {
    #gam_fit <- mgcv::gam(Y ~ s(exposure) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
    gam_fit <- mgcv::gam(Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                        data = data_example, family = family)
  } else {
    gam_fit <- gam::gam(Y ~ s(exposure, bs = 'cr', k = 4), data = data_example, family = family)
  }
  
  # Now fit the causal model with weight by GPS
    
    # Create dataset for GPS estimation
    data_gps <- dplyr::select(data_example, -Y)
  
    # Generate weights using entropy balancing weights 
    W1 <- weightit(exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_gps, 
                   method = "ebal", stabilize = T)
    data_gps$ent <- W1$weight
    
    # Truncate the 99th percentile
    upper_bound = quantile(data_gps$ent, 0.995)
    data_gps[ent > upper_bound, ent := upper_bound]
    
    # Generate weights balancing on second moment 
    W2 <- weightit(exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_gps, 
                   method = "ebal", stabilize = T, moments = 2)
    data_gps$ent2 <- W2$weight
    
    # Truncate the 99th percentile
    upper_bound = quantile(data_gps$ent2, 0.995)
    data_gps[ent2 > upper_bound, ent2 := upper_bound]
    
    # Generate weights based on "energy" balancing
    W3 <- weightit(exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_gps, 
                   method = "energy")
    data_gps$ent3 <- W3$weight
    
    # Merge back on entropy weights to original dataset
    data_ent <- 
      data_example %>% 
      left_join(data_gps, by = c("exposure", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6"))
    
    # Create correlation table
    correlation_table <- 
      rbindlist(lapply(c("ent", "ent2", "ent3"), function(entropy_type) {
        # Grab weight of interest
        correct_weight <- pull(data_ent, !!entropy_type)
        post_cor <- cov.wt(data_ent[, c("exposure", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6")], 
                           wt = correct_weight, cor = T)$cor
        post_cor <- abs(post_cor[-1, 1])
        
        pre_cor <- cov.wt(data_ent[, c("exposure", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6")], cor = T)$cor
        pre_cor <- abs(pre_cor[-1, 1])
        
        correlation_table <- data.table(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"), 
                                        pre_cor = pre_cor, 
                                        post_cor = post_cor, 
                                        method = entropy_type)
      }))
   
    
    
    # Now fit using CausalGPS package --------------------------
    if (causal_gps) {
      # First fit using default settings for our simulation
      pseudo_pop_default <- generate_pseudo_pop(data_ent$Y,
                                                data_ent$exposure,
                                                data.frame(data_ent[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transform = TRUE,
                                                transformers = list("pow2", "pow3", "abs", "scale"),
                                                trim_quantiles = c(0, 1),
                                                optimized_compile = TRUE,
                                                sl_lib = c("m_xgboost"),
                                                covar_bl_method = "absolute",
                                                covar_bl_trs = 0.1,
                                                covar_bl_trs_type = "mean",
                                                max_attempt = 1,
                                                matching_fun = "matching_l1",
                                                delta_n = 1.0,
                                                scale = 1,
                                                nthread = 2)
      
      # Now fit semi-parametric here 
      causal_gps_default <- 
        mgcv::gam(formula = Y ~ s(w, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                  family = "gaussian",
                  data = data.frame(pseudo_pop_default$pseudo_pop),
                  weights = counter_weight)
      
      # Now fit most optimized causal pathway
      # Set grid to tune over 
      #nrounds = 300
      tune_grid <- expand.grid(
        nrounds = c(100),
        eta = c(0.05, 0.1, 0.2, 0.3),
        max_depth = c(2, 3),
        delta = seq(1, 3, by = 0.1)
      )
      
      # Make list to pass to run in parallel using mclapply
      tune_grid_list <- as.list(as.data.frame(t(tune_grid)))
      
      # Wrapper function for running causalGPS
      wrapper_func <- function(tune_param){
        pseudo_pop_tune <- generate_pseudo_pop(data_ent$Y,
                                               data_ent$exposure,
                                               data.frame(data_ent[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                               ci_appr = "matching",
                                               pred_model = "sl",
                                               gps_model = "parametric",
                                               use_cov_transform = FALSE,
                                               transformers = list("pow2", "pow3", "abs", "scale"),
                                               trim_quantiles = c(0.01, 0.99),
                                               optimized_compile = T,
                                               sl_lib = c("m_xgboost"),
                                               params = list(xgb_rounds = tune_param[[1]],
                                                             xgb_eta = tune_param[[2]],
                                                             xgb_max_depth = tune_param[[3]]),
                                               covar_bl_method = "absolute",
                                               covar_bl_trs = 0.1,
                                               covar_bl_trs_type = "mean",
                                               max_attempt = 1,
                                               matching_fun = "matching_l1",
                                               delta_n = tune_param[[4]],
                                               scale = 1,
                                               nthread = 1)
        
        
        matched_pop_tune <- pseudo_pop_tune$pseudo_pop
        
        # Truncate upper 1%
        # matched_pop[counter_weight > quantile(matched_pop$counter_weight, 0.99),
        #             counter_weight := quantile(matched_pop$counter_weight, 0.99)]
        
        # Now generate covariate balance tab
        balance_table <- 
          bal.tab(w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                  data = matched_pop_tune,
                  weights = matched_pop_tune$counter_weight,
                  method = "weighting",
                  stats = c("cor"),
                  un = T,
                  continuous = "std",
                  s.d.denom = "weighted",
                  abs = T,
                  thresholds = c(cor = .1), poly = 1)
        
        mean_post_cor <- mean(balance_table$Balance$Corr.Adj)
        
        results <- 
          data.frame(nrounds = tune_param[[1]], 
                     eta = tune_param[[2]],
                     max_depth = tune_param[[3]],
                     delta = tune_param[[4]],
                     scale = 1,
                     mean_corr = mean_post_cor,
                     max_corr = max(balance_table$Balance$Corr.Adj))
        return(list(results, matched_pop_tune))
        
      }
      
      # Now iterate through simulation function, bind together results
      pseudo_pop_list <- mclapply(tune_grid_list, mc.cores = 6, wrapper_func)
      corr_search <- do.call("rbind", (lapply(pseudo_pop_list, function(x) x[[1]])))
      min_corr = corr_search %>% filter(mean_corr == min(mean_corr)) %>% slice_head()
      
      # Extract minimum as your result
      pseudo_pop_tuned <- wrapper_func(min_corr[1:4])[[2]]
      
      # Check that correlation
      post_cor <- cov.wt(pseudo_pop_tuned[, c("w", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6")], wt = (pseudo_pop_tuned$counter_weight), cor = T)$cor
      post_cor <- abs(post_cor[-1, 1])
      
      # Now fit semi-parametric here 
      causal_gps_tuned <- 
        mgcv::gam(formula = Y ~ s(w, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                  family = "gaussian",
                  data = data.frame(pseudo_pop_tuned),
                  weights = counter_weight)
    }
    
    # Dropping these models now 
    if (adjust_confounder) {
      entropy_lm <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                           data = data_ent, family = family, weights = ent)
      entropy_gam <- mgcv::gam(Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_ent, family = family, weights = ent)
      
      # Now fit second moment entropy weighting
      entropy_lm2 <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                           data = data_ent, family = family, weights = ent2)
      entropy_gam2 <- mgcv::gam(Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                                data = data_ent, family = family, weights = ent2)
      
      # Now fit "energy" 
      entropy_lm3 <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                         data = data_ent, family = family, weights = ent3)
      entropy_gam3 <- mgcv::gam(Y ~ s(exposure, bs = 'cr', k = 4) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, 
                                data = data_ent, family = family, weights = ent3)
      
    } else {
      propensity_lm <- glm(Y ~ exposure, data = data_ent, family = family, weights = ent)
      propensity_gam <- gam::gam(Y ~ s(exposure, df = 3), data = data_ent, family = family, weights = ent)
    }
    
    # Now fit threshold model 
    if (adjust_confounder) {
      change_model <- chngptm(Y ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, ~ exposure, family = "gaussian", type = "segmented", var.type = "bootstrap", save.boot = T, data = data_ent)
      change_model_ent <- chngptm(Y ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, ~ exposure, family = "gaussian", 
                                  type = "segmented", var.type = "bootstrap", save.boot = T, data = data_ent, weights = data_ent$ent)
      
    } else {
      # fit unadjusted model
      change_model <- chngptm(Y ~ 1, ~ exposure, family = "gaussian", type = "segmented", var.type = "bootstrap", save.boot = T, data = data_example)
    }
    
    # Lets extract what the proposed change model threshold is
    change_assesed = change_model$chngpt
    change_lower = summary(change_model)$chngpt[["(lower"]]
    change_upper = summary(change_model)$chngpt[["upper)"]]
    
    # Do the same change point assessed using entropy weights 
    change_assesed_ent = change_model_ent$chngpt
    change_lower_ent = summary(change_model_ent)$chngpt[["(lower"]]
    change_upper_ent = summary(change_model_ent)$chngpt[["upper)"]]
  
  # Now many points to evaluate the integral
  discrete_points = 100
  
  # Predict the RC curve for every point
  ## Fix this for nonlinear first
  data_prediction <- 
    rbindlist(lapply(seq(0, 20, length.out = discrete_points), function(pot_exp) {
      potential_data <- 
        dplyr::select(data_example, cf1, cf2, cf3, cf4, cf5, cf6) %>% 
        mutate(exposure = pot_exp)
      
      # Add on exposure effect depending on relationship
      if (exposure_relationship == "linear") {
        potential_data$transformed_exp <- potential_data$exposure
        Y = potential_data$transformed_exp
      } else if (exposure_relationship == "sublinear") {
        potential_data$transformed_exp = 5 * log(potential_data$exposure + 1)
        Y = potential_data$transformed_exp
      } else if (exposure_relationship == "threshold") {
        # Set threshold at 5
        thresh_exp <- potential_data$exposure 
        thresh_exp[thresh_exp <= 5] <- 0
        thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
        potential_data$transformed_exp <- 1.5 * thresh_exp
        Y = potential_data$transformed_exp
      }
      
      # Add on effect of the confounders
      if (outcome_relationship == "linear") {
        Y = Y + as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(dplyr::select(potential_data, cf1, cf2, cf3, cf4, cf5, cf6)))
      } else if (outcome_relationship == "interaction") {
        # problem with the threshold exposure scenario happening 
        Y = Y + as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(dplyr::select(potential_data, cf1, cf2, cf3, cf4, cf5, cf6)) - 
                             potential_data$transformed_exp*(-0.1 * potential_data$cf1 + 0.1 * potential_data$cf3^2 + 0.1*potential_data$cf4 + 0.1*potential_data$cf5))
      } else {
        stop("Outcome relationship for now is only between linear and interaction")
      }
      
      # Fit each model and take mean for ERC
      potential_outcome <- 
        potential_data %>% 
        mutate(linear_model = predict(linear_fit, potential_data, type = "response"),
               gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
               ent_linear = predict(entropy_lm, newdata = potential_data, type = "response"),
               ent_gam = predict(entropy_gam, newdata = potential_data, type = "response"),
               ent_linear2 = predict(entropy_lm2, newdata = potential_data, type = "response"), # Adjust for second moment
               ent_gam2 = predict(entropy_gam2, newdata = potential_data, type = "response"),
               ent_linear3 = predict(entropy_lm3, newdata = potential_data, type = "response"), # Energy balancing
               ent_gam3 = predict(entropy_gam3, newdata = potential_data, type = "response"),
               change_model = predict(change_model, newdata = potential_data, type = "response"),
               change_ent = predict(change_model_ent, newdata = potential_data, type = "response"),
               causal_gps_tuned = potential_data %>% rename(w = exposure) %>% predict(causal_gps_tuned, newdata = ., type = "response"),
               true_fit = Y) %>% 
        dplyr::select(exposure, linear_model, gam_model, ent_linear, ent_gam, ent_linear2, ent_gam2,
                      ent_linear3, ent_gam3, change_model, change_ent, causal_gps_tuned, true_fit) %>% 
        summarize_all(mean) %>% 
        data.table()
      return(potential_outcome)
    }))
  
  model_types <- c("linear_model", "gam_model", "ent_linear", "ent_gam", "ent_linear2", "ent_gam2",
                   "ent_linear3", "ent_gam3", "causal_gps_tuned", "change_model", "change_ent")
  
  # For now evaluting curve from 
  trim_upper <- 20
  data_prediction <- data_prediction[exposure <= trim_upper]
  
  # Diagnostic plot of model fit
  plot_fits <- 
    data_prediction %>% 
    pivot_longer(c(all_of(model_types), "true_fit"), names_to = "model", values_to = "prediction") %>%
    arrange(exposure) %>%
    ggplot(aes(x = exposure, y = prediction, color = model, linetype = model)) +
    geom_line() +
    labs(x = "Exposure concentration", y = "Relative risk of death")
  
  # Add specified change point here for both unweighted and entropy weighted example
  data_prediction$change_point = change_assesed
  data_prediction$change_lower = change_lower
  data_prediction$change_upper = change_upper
  
  data_prediction$change_point_ent = change_assesed_ent
  data_prediction$change_lower_ent = change_lower_ent
  data_prediction$change_upper_ent = change_upper_ent
  
  
  data_metrics <- 
    data_prediction %>% 
    tidyr::gather(all_of(model_types), key = "model", value = "prediction") %>% 
    data.table()

  # Changing metric to calculate the bias, and also the absolute bias. Return it all for now
  data_metrics <- 
    data_metrics[, .(bias = prediction - true_fit,
                     mse = (prediction - true_fit) ^ 2), by = .(model, exposure)]
  
  return_columns <- c("exposure", model_types, "true_fit", "change_point", "change_lower", "change_upper",
                      "change_point_ent", "change_lower_ent", "change_upper_ent")
  return(list(metrics = data_metrics, predictions = data_prediction[, ..return_columns], cor_table = correlation_table, 
              pseudo_pop_tuned = pseudo_pop_tuned))
}

