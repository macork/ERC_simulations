# File that contains all the functions I am using for the simulation 

data_generate_a <-
  function(sample_size = 1000, exposure = NA, confounders = NA, exposure_relationship = NA, outcome_relationship = NA, family = NA, outcome_sd = 10) {
    
    # Abbreviate for clarity 
    cf <- confounders
    
    # For now mostly using gaussian family, though poisson has been used in the past as well
    if (family == "gaussian") {
      
      # Add on exposure effect depending on relationship
      if (exposure_relationship == "linear") {
        transformed_exp <- exposure
        Y = transformed_exp
      } else if (exposure_relationship == "sublinear") {
        transformed_exp = log10(exposure + 1)
        Y = 8 * transformed_exp
      } else if (exposure_relationship == "threshold") {
        # Set threshold at 5
        transformed_exp <- exposure 
        transformed_exp[transformed_exp <= 5] <- 0
        transformed_exp[transformed_exp > 5] <- transformed_exp[transformed_exp > 5] - 5
        Y = Y + transformed_exp
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

# Define function for simulations
# adjust_confounder says if we should adjust for confounders at all in our model
metrics_from_data <- function(exposure = NA, exposure_relationship = "linear", outcome_relationship = "linear", sample_size = 1000, confounders = NA,
                              family = "gaussian", eschif_draws = NULL, adjust_confounder = T, causal_gps = F) {
  
  
  # exposure_relationship = "linear"
  # outcome_relationship = "interaction"
  #  family = "gaussian"
  #  adjust_confounder = T
  
  # Simulate data
  data_example <- data_generate_a(sample_size = sample_size, exposure = exposure, confounders = confounders, 
                                  exposure_relationship = exposure_relationship, outcome_relationship = outcome_relationship,
                                  family = family)
  
  # Now fit a linear model
  if (adjust_confounder == T) {
    linear_fit <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
  } else {
    linear_fit <- glm(Y ~ exposure, data = data_example, family = family)
  }
  
  # Now fit a GAM model 
  if (adjust_confounder == T) {
    #gam_fit <- mgcv::gam(Y ~ s(exposure) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
    gam_fit <- gam::gam(Y ~ s(exposure, df = 3) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
  } else {
    gam_fit <- gam::gam(Y ~ s(exposure, df = 3), data = data_example, family = family)
  }
  
  # Now fit the causal model with weight by GPS
    
    # Create dataset for GPS estimation
    data_gps <- dplyr::select(data_example, -Y)
    
    # # Linear GPS model
    # GPS_mod <- lm(exposure ~ ., data = data_gps)
    # # Create stabilized ipw weight from linear GPS
    # ipw <- dnorm(data_gps$exposure, mean = mean(data_gps$exposure), sd = sd(data_gps$exposure)) / dnorm(data_gps$exposure, predict(GPS_mod, data_gps), sqrt(mean(resid(GPS_mod)^2)))
    # data_gps$IPW <- ipw
  
    # Generate weights using entropy balancing weights 
    W1 <- weightit(exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_gps, method = "ebal", stabilize = T)
    data_gps$IPW <- W1$weights
    
    
    # Train Xgboost for predicting exposure (not done for now)
    # if (causal_gps) {
    #   # The following lines are using caret but I am going down to two options: "default" and "tuned" which is the lowest average correlation
    #   tune_control <-
    #     caret::trainControl(method = "cv",
    #                         number = 5,
    #                         verboseIter = FALSE,
    #                         allowParallel = TRUE)
    #   nrounds <- 700
    #   tune_grid <- expand.grid(nrounds = seq(from = 400, to = nrounds, by = 50),
    #                            eta = c(0.025, 0.05, 0.1),
    #                            max_depth = c(2, 3),
    #                            gamma = 0,
    #                            colsample_bytree = 1,
    #                            min_child_weight = c(1),
    #                            subsample = 1)
    #   xgbt_model <-
    #     train(exposure ~ ., data = data_gps, method = "xgbTree",
    #           trControl = tune_control, tuneGrid = tune_grid, verbosity = 0)
    # 
    #   xgbt_model$bestTune
    #   
    #   predVals <- extractPrediction(list(xgbt_model))
    #   plotObsVsPred(predVals)
    #   
    #   ipw <- dnorm(data_gps$exposure, mean(data_gps$exposure), sd(data_gps$exposure)) / dnorm(data_gps$exposure, predVals$pred, sqrt(mean(resid(xgbt_model)^2)))
    #   data_gps$IPW <- ipw
    # }

    
    # Truncate the 99th percentile (could also truncate to 10 as Xiao did in his paper)
    upper_bound = quantile(data_gps$IPW, 0.99)
    # lower_bound = quantile(data_gps$IPW, 0.01)
    
    # Truncating above 10 now (following Xiao advice)
    data_gps[IPW > 10, IPW := 10]
    
    data_ipw <- data_example %>% left_join(data_gps, by = c("exposure", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6"))
    #data_ipw = filter(data_example, between(IPW, lower_bound, upper_bound))
    
    # Make correlation table.... might be difficult 
    
    post_cor <- cov.wt(data_ipw[, c("exposure", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6")], wt = (data_ipw$IPW), cor = T)$cor
    post_cor <- abs(post_cor[-1, 1])
    
    pre_cor <- cov.wt(data_ipw[, c("exposure", "cf1", "cf2", "cf3", "cf4", "cf5", "cf6")], cor = T)$cor
    pre_cor <- abs(pre_cor[-1, 1])
    
    correlation_table <- data.table(covariate = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6"), pre_cor = pre_cor, post_cor = post_cor)
    
    # correlation_table <- 
    #     rbindlist(lapply(names(data.frame(confounders)), function(cov){
    #       pre_cor <- abs(cor(data_ipw$exposure, data_ipw[[cov]]))
    #       post_cor <- abs(wtd.cor(data_ipw$exposure, data_ipw[[cov]], data_ipw$IPW))
    #       data.table(covariate = cov, pre_cor = pre_cor, post_cor = post_cor[1])
    #     }))
    
    
    # Now fit using CausalGPS package --------------------------
    if (causal_gps) {
      # First fit using default settings for our simulation
      pseudo_pop_default <- generate_pseudo_pop(data_ipw$Y,
                                        data_ipw$exposure,
                                        data.frame(data_ipw[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                        ci_appr = "matching",
                                        pred_model = "sl",
                                        gps_model = "parametric",
                                        use_cov_transform = TRUE,
                                        transformers = list("pow2", "pow3", "abs", "scale"),
                                        trim_quantiles = c(0.025,0.975),
                                        optimized_compile = TRUE,
                                        sl_lib = c("m_xgboost"),
                                        covar_bl_method = "absolute",
                                        covar_bl_trs = 0.1,
                                        covar_bl_trs_type = "mean",
                                        max_attempt = 20,
                                        matching_fun = "matching_l1",
                                        delta_n = 1.0,
                                        scale = 1,
                                        nthread = 2)
      
      # Now fit semi-parametric here 
      causal_gps_default <- gam::gam(formula = Y ~ s(w, df = 3),
                             family = "gaussian", data = data.frame(pseudo_pop_default$pseudo_pop), 
                             weights = counter)
      
      # Now fit most optimized causal pathway
      # Set grid to tune over 
      nrounds = 300
      tune_grid <- expand.grid(
        nrounds = nrounds,
        eta = c(0.025, 0.05, 0.1),
        max_depth = c(2, 3),
        delta = seq(0.2, 2, by = 0.2)
      )
      
      # Make list to pass to parLapply
      tune_grid_list <- as.list(as.data.frame(t(tune_grid)))
      
      # Wrapper function for running causalGPS
      wrapper_func <- function(tune_grid){
        pseudo_pop <- generate_pseudo_pop(data_ipw$Y,
                                          data_ipw$exposure,
                                          data.frame(data_ipw[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]),
                                          ci_appr = "matching",
                                          pred_model = "sl",
                                          gps_model = "parametric",
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
                                          nthread = 2)
        
        
        results <- 
          data.frame(nrounds = tune_grid[[1]], 
                     eta = tune_grid[[2]],
                     max_depth = tune_grid[[3]],
                     delta = tune_grid[[4]],
                     scale = 1,
                     mean_corr = pseudo_pop$adjusted_corr_results$mean_absolute_corr,
                     max_corr = pseudo_pop$adjusted_corr_results$maximal_absolute_corr)
        return(list(results, pseudo_pop))
        
      }
      
      # Run these wrapper functions in parallel
      # Consider running optimized vs non optimized version through simulations? Could try to get that working today 
      
      # currently not running parallel so commenting this out 
      # cl <- parallel::makeCluster(10, type="PSOCK")
      # parallel::clusterExport(cl=cl,
      #                         varlist = c("generate_pseudo_pop",
      #                                     "wrapper_func",
      #                                     "data_ipw"
      #                         ),
      #                         envir=environment())
      # 
      # pseudo_pop_list  <- parallel::parLapply(cl,tune_grid_list, wrapper_func)
      # parallel::stopCluster(cl)
      
      pseudo_pop_list <- lapply(tune_grid_list, wrapper_func)
      
      # If you don't want to run in parallel, use this 
      #pseudo_pop_list <- mclapply(tune_grid_list, mc.cores = 4, wrapper_func)
      corr_search <- do.call("rbind", (lapply(pseudo_pop_list, function(x) x[[1]])))
      
      # Extract minimum as your result
      pseudo_pop_tuned <- pseudo_pop_list[[which.min(corr_search$mean_corr)]][[2]]
      
      # Now fit semi-parametric here 
      causal_gps_tuned <- gam::gam(formula = Y ~ s(w, df = 3), family = "gaussian", data = data.frame(pseudo_pop_tuned$pseudo_pop), weights = counter)
      
    }
    
    # Dropping these models now 
    if (adjust_confounder) {
      propensity_lm <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_ipw, family = family, weights = IPW)
      propensity_gam <- gam::gam(Y ~ s(exposure, df = 3) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_ipw, family = family, weights = IPW)
    } else {
      propensity_lm <- glm(Y ~ exposure, data = data_ipw, family = family, weights = IPW)
      propensity_gam <- gam::gam(Y ~ s(exposure, df = 3), data = data_ipw, family = family, weights = IPW)
    }
    
    # Propensity gam being supressed for now (not sure if this is correct formulation)
    #propensity_gam <- gam::gam(Y ~ s(exposure, df = 3) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_ipw, family = family, weights = IPW)
    #propensity_gam_nocf <- mgcv::gam(Y ~ s(exposure), data = data_ipw, family = family, weights = IPW)
    
    # Now fit threshold model 
    if (adjust_confounder) {
      change_model <- chngptm(Y ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, ~ exposure, family = "gaussian", type = "segmented", var.type = "bootstrap", save.boot = T, data = data_example)
    } else {
      # fit unadjusted model
      change_model <- chngptm(Y ~ 1, ~ exposure, family = "gaussian", type = "segmented", var.type = "bootstrap", save.boot = T, data = data_example)
    }
    
    # Suppress propensity change point model for now, seemingly not working
    #propensity_change <- chngptm(Y ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, ~ exposure, family = "gaussian", type = "segmented", var.type = "bootstrap", save.boot = T, data = data_ipw, weights = data_ipw$IPW)
    
    # Lets extract what the proposed change model threshold is
    change_assesed = change_model$chngpt
    change_lower = summary(change_model)$chngpt[["(lower"]]
    change_upper = summary(change_model)$chngpt[["upper)"]]
  
  # Now many points to evaluate the integral
  discrete_points = 1000
  
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
        potential_data$transformed_exp = log10(potential_data$exposure + 1)
        Y = 8 * potential_data$transformed_exp
      } else if (exposure_relationship == "threshold") {
        # Set threshold at 5
        thresh_exp <- potential_data$exposure 
        thresh_exp[thresh_exp <= 5] <- 0
        thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
        potential_data$transformed_exp <- thresh_exp
        Y = Y + thresh_exp
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
      
      # # Add on exposure effect
      # if (exposure_relationship == "linear") {
      #   Y = Y + potential_data$exposure
      # } else if (exposure_relationship == "sublinear") {
      #   Y = Y + 8 * log10(potential_data$exposure + 1)
      # } else if (exposure_relationship == "threshold") {
      #   # Set threshold at 5
      #   thresh_exp <- potential_data$exposure 
      #   thresh_exp[thresh_exp <= 5] <- 0
      #   thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
      #   Y = Y + thresh_exp
      # }
      
      # Fit each model and take mean for ERC
      potential_outcome <- 
        potential_data %>% 
        mutate(linear_model = predict(linear_fit, potential_data, type = "response"),
               gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
               linear_gps = predict(propensity_lm, newdata = potential_data, type = "response"),
               gam_gps = predict(propensity_gam, newdata = potential_data, type = "response"),
               change_model = predict(change_model, newdata = potential_data, type = "response"),
               causal_gps_default = potential_data %>% rename(w = exposure) %>% predict(causal_gps_default, newdata = ., type = "response"),
               causal_gps_tuned = potential_data %>% rename(w = exposure) %>% predict(causal_gps_tuned, newdata = ., type = "response"),
               true_fit = Y) %>% 
        dplyr::select(exposure, linear_model, gam_model, linear_gps, gam_gps, change_model, causal_gps_default, causal_gps_tuned, true_fit) %>% 
        summarize_all(mean) %>% 
        data.table()
      return(potential_outcome)
    }))
  
  
  # Here add manual code to say that for the changepoint you can just add the linear prediction after it thinks it found a change
  # threshold = summary(change_model)$chngpt[["est"]]
  # coef_exposure_pre = coef(change_model)[["exposure"]]
  # coef_exposure_post = coef(change_model)[["(exposure-chngpt)+"]]
  # intercept = coef(change_model)[["(Intercept)"]]
  # potential_outcome$change_model = ifelse(potential_data$exposure < threshold, intercept + coef_exposure_pre*(potential_data$exposure), intercept + coef_exposure_pre*(potential_data$exposure) + coef_exposure_post*(potential_data$exposure - threshold))
  
  # Add changepoint for propensity model, suppressed for now since didn't seem to improve things 
  # threshold = summary(propensity_change)$chngpt[["est"]]
  # coef_exposure_pre = coef(propensity_change)[["exposure"]]
  # coef_exposure_post = coef(propensity_change)[["(exposure-chngpt)+"]]
  # intercept = coef(propensity_change)[["(Intercept)"]]
  # data_prediction$propensity_change = ifelse(data_prediction$exposure < threshold, intercept + coef_exposure_pre*(data_prediction$exposure), intercept + coef_exposure_pre*(data_prediction$exposure) + coef_exposure_post*(data_prediction$exposure - threshold))
  
  # Calculate trimmed reference exposure value (set as minimum so that the curve starts at 0)
  trimmed_reference <- sort(data_prediction$exposure)[1]
  
  # Remove influence of the intercept to just compare risk difference or relative risk 
  if (family == "gaussian") {
      # add GAM GPS here if using that model 
    if (causal_gps) {
      model_types <- c("linear_model", "gam_model", "linear_gps", "gam_gps", "causal_gps_default", "causal_gps_tuned", "change_model")
      data_prediction <- 
        data_prediction %>% 
        filter(exposure >= trimmed_reference) %>% 
        mutate(linear_model = linear_model - data_prediction[exposure == trimmed_reference, linear_model],
               gam_model = gam_model - data_prediction[exposure == trimmed_reference, gam_model],
               linear_gps = linear_gps - data_prediction[exposure == trimmed_reference, linear_gps],
               gam_gps = gam_gps - data_prediction[exposure == trimmed_reference, gam_gps],
               change_model = change_model - data_prediction[exposure == trimmed_reference, change_model],
               causal_gps_default = causal_gps_default - data_prediction[exposure == trimmed_reference, causal_gps_default],
               causal_gps_tuned = causal_gps_tuned - data_prediction[exposure == trimmed_reference, causal_gps_tuned],
               #propensity_change = propensity_change - data_prediction[exposure == trimmed_reference, propensity_change],
               true_fit = true_fit - data_prediction[exposure == trimmed_reference, true_fit])
    } else {
      model_types <- c("linear_model", "gam_model", "linear_gps", "gam_gps", "change_model")
      data_prediction <- 
        data_prediction %>% 
        filter(exposure >= trimmed_reference) %>% 
        mutate(linear_model = linear_model - data_prediction[exposure == trimmed_reference, linear_model],
               gam_model = gam_model - data_prediction[exposure == trimmed_reference, gam_model],
               linear_gps = linear_gps - data_prediction[exposure == trimmed_reference, linear_gps],
               gam_gps = gam_gps - data_prediction[exposure == trimmed_reference, gam_gps],
               change_model = change_model - data_prediction[exposure == trimmed_reference, change_model],
               #propensity_change = propensity_change - data_prediction[exposure == trimmed_reference, propensity_change],
               true_fit = true_fit - data_prediction[exposure == trimmed_reference, true_fit])
    }
  } else if (family == "poisson") {
    model_types <- c("linear_model", "gam_model")
    data_prediction <- 
      data_prediction %>% 
      mutate(linear_model = linear_model / data_prediction$linear_model[1],
             gam_model = gam_model / data_prediction$gam_model[1],
             true_fit = true_fit / data_prediction$true_fit[1])
    
    # Currently not using eSCHIF 
    if ("eschif" %in% model_types) {
    # Now fit eSCHIF ---------------------------------
    range <- max(exposure) - min(exposure)
    alpha = seq(1, range, by = 2)
    mu = seq(0, range, by = 2)
    tau = c(0.1, 0.2, 0.4, 0.8, 1)
    thres = seq(0, range, 0.5)
    
    # Get best fit for eSCHIF
    # start.time <- Sys.time()
    if (is.null(eschif_draws)) {
      y_eschif <- log(data_prediction$gam_model)
      eschif_fit <-
        rbindlist(lapply(alpha, function(a) {
          rbindlist(lapply(mu, function(m) {
            rbindlist(lapply(tau, function(t) {
              rbindlist(lapply(thres, function(th) {
                z = ((exposure - th) + abs(exposure - th)) / 2
                diff = log(z / a + 1) / (1 + exp(-(z - m) / (t * range)))
                fit = lm(y_eschif ~ diff - 1)
                data.table(alpha = a, mu = m, tau = t, thres = th, aic = AIC(fit), theta = coef(fit))
              }))
            }))
          }))[aic == min(aic)]
        }))[aic == min(aic)]
      # end.time <- Sys.time()
      # time.taken <- end.time - start.time
      # time.taken
      
      z = ((exposure - eschif_fit$thres) + abs(exposure - eschif_fit$thres)) / 2
      eschif_pred = exp(eschif_fit$theta * log(z / eschif_fit$alpha + 1) / (1 + exp(-(z - eschif_fit$mu) / (eschif_fit$tau * range))))
      data_prediction$eschif <- eschif_pred
    } else {
      # Run through eSCHIF as many times as specified 
      stop("Not done with this part of function yet")
      
    }
    }
  }
  
  # For now try trimming the 10% boundary at the top (wonder what effect this has, just setting to 15 for now)
  trim_upper <- 15
  data_prediction <- data_prediction[exposure <= trim_upper]
  
  
  # Diagnostic plot of model fit
  plot_fits <- 
    data_prediction %>% 
    pivot_longer(c(all_of(model_types), "true_fit"), names_to = "model", values_to = "prediction") %>%
    arrange(exposure) %>%
    ggplot(aes(x = exposure, y = prediction, color = model, linetype = model)) +
    geom_line() +
    labs(x = "Exposure concentration", y = "Relative risk of death")
  
  # Add specified change point here
  data_prediction$change_point = change_assesed
  data_prediction$change_lower = change_lower
  data_prediction$change_upper = change_upper
  
  
  data_metrics <- 
    data_prediction %>% 
    tidyr::gather(model_types, key = "model", value = "prediction") %>% 
    data.table()
  
  data_metrics <- 
    data_metrics[, .(bias = prediction - true_fit,
                     mse = (prediction - true_fit) ^ 2), by = .(model)]
  
  data_metrics <- data_metrics[, .(bias = mean(bias), mse = mean(mse)), by = .(model)]
  return_columns <- c("exposure", model_types, "true_fit", "change_point", "change_lower", "change_upper")
  return(list(metrics = data_metrics, predictions = data_prediction[, ..return_columns], cor_table = correlation_table, 
              pseudo_pop_default = pseudo_pop_default, pseudo_pop_tuned = pseudo_pop_tuned))
}


# Function for evaluating eSCHIF 
## Currently that will take forever to run because eSCHIF runs 1,000 versions of its function anyway, so thinking if that is worth it. 




