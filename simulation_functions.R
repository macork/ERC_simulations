# File that contains all the functions I am using for the simulation 

data_generate_a <-
  function(sample_size = 1000, exposure = NA, confounders = NA, relationship = NA, family = NA, outcome_sd = 10, confounder_mult) {
    
    # Abbreviate for clarity 
    cf <- confounders
    
    if (relationship == "linear") {
      if (family == "gaussian") {
        Y = as.numeric(20 - confounder_mult * c(2, 2, 3, -1, -2, -2) %*% t(cf) + 0.1 * exposure + rnorm(sample_size, mean = 0, sd = outcome_sd))
      } else if (family == "poisson") {
        # Poisson model
        lambdas = exp(2 +  as.vector(matrix(c(0.2, 0.2, 0.3, -0.1, -0.2, 0.2), nrow = 1) %*% t(cf)) + 0.1 * exposure)
        Y = rpois(n = 1000, lambda = lambdas)
      }
      
      #Y_true = as.numeric(20 - c(2, 2, 3, -1) %*% t(cf) - 2 * cf5 - 2 * cf6 + 0.1 * exposure)
    } else if (relationship == "sublinear") {
      #Y = as.numeric(20 - c(2, 2, 3, -1) %*% t(cf) - 2 * cf5 - 2 * cf6 + 5 * sqrt(exposure) + rnorm(size, mean = 0, sd = outcome_sd))
      if (family == "gaussian") {
        Y = as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(cf) + 8 * log10(exposure + 1) + rnorm(sample_size, mean = 0, sd = outcome_sd))
        #Y_true = as.numeric(20 - c(2, 2, 3, -1) %*% t(cf) - 2 * cf5 - 2 * cf6 + log(exposure))
      } else if (family == "poisson") {
        lambdas = exp(as.numeric(unlist(2 + 0.2 * cf[, 1] + 0.2 * cf[, 2] + 0.3 * cf[, 3] -0.1 * cf[, 4] - 0.2 * cf[, 5] + 0.2 * cf[, 6] + 1.5 * log10(exposure + 1))))
        Y = rpois(n = 1000, lambda = lambdas)
      }
    
    } else if (relationship == "threshold") {
      if (family == "gausian") {
        thresh_exp <- exposure 
        thresh_exp[thresh_exp <= 5] <- 0
        thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
        Y = as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(cf) + 1 * thresh_exp + rnorm(sample_size, mean = 0, sd = outcome_sd))
      } else if (family == "poisson") {
        thresh_exp <- exposure 
        thresh_exp[thresh_exp <= 5] <- 0
        thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
        lambdas = exp(as.numeric(unlist(2 + 0.2 * cf[, 1] + 0.2 * cf[, 2] + 0.3 * cf[, 3] -0.1 * cf[, 4] - 0.2 * cf[, 5] + 0.2 * cf[, 6] + 0.1 * thresh_exp)))
        Y = rpois(n = 1000, lambda = lambdas)
      }
    } else {
      stop("Must either be linear, sublinear, or threshold as of now")
    }
    
    
    # Save into one simulated dataset
    simulated_data <- data.table(cbind(Y, exposure, cf))
    colnames(simulated_data)[3:8] <- c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")
    return(simulated_data)
  }

# Define function for simulations
# adjust_confounder says if we should adjust for confounders at all in our model
metrics_from_data <- function(just_plots = F, exposure = NA, relationship = "linear", sample_size = 1000,confounders = NA,
                              family = "gaussian", eschif_draws = NULL, adjust_confounder = T, confounder_mult = 1) {
  
  # Fit data generating mechanism
  data_example <- data_generate_a(sample_size = sample_size, exposure = exposure, confounders = confounders, 
                                  relationship = relationship, family = family, confounder_mult = confounder_mult)
  
  # Now fit a linear model
  if (adjust_confounder == T) {
    linear_fit <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
  } else {
    linear_fit <- glm(Y ~ exposure, data = data_example, family = family)
  }
  
  # Now fit a GAM model 
  if (adjust_confounder == T) {
    gam_fit <- mgcv::gam(Y ~ s(exposure) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
  } else {
    gam_fit <- mgcv::gam(Y ~ s(exposure), data = data_example, family = family)
  }
  
  # Now fit the causal model with weight by GPS
  if (adjust_confounder == T) {
    
    # Create dataset for GPS
    data_gps <- data.table(confounders, exposure)
    
    # # Linear GPS model, skipping for now
    # GPS_mod <- lm(exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6 + cf1^2 + cf2^2 + cf3^2 + cf4^2 + cf5^2 + cf6^2, data = data_gps)
    # mean_predict <- predict(GPS_mod3, data_gps)
    # 
    # GPS_mod <- CBPS(exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_gps)

    # Different propensity score measures
    W1 <- weightit(exposure ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6 + cf1^2 + cf2^2 + cf3^2 + cf4^2 + cf5^2 + cf6^2, data = data_gps, method = "npcbps", estimand = "ATT", over = F)
    data_example$IPW <- W1$weights
    
    
    
    
    # Or nonlinear GPS model
    # sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.glmnet", "SL.ranger")
    # GPS_mod <- suppressWarnings(SuperLearner(Y = exposure, X = data.frame(confounders), SL.library = sl.lib, family = gaussian()))
    # mean_predict <- c(GPS_mod$SL.predict)
    # mod_sd <- sd(exposure - mean_predict)
    # feature_names <- GPS_mod$feature_names
    # GPS <- dnorm(exposure, mean = mean_predict, sd = mod_sd) # probability of exposure given confounders
    # Nm <- dnorm(exposure, mean = mean(exposure, na.rm = TRUE), sd = sd(exposure, na.rm = TRUE)) # General probability of exposure
    # data_example$IPW <- Nm / GPS
    

    # cv_sl = CV.SuperLearner(Y = exposure, X = data.frame(confounders), family = gaussian(),
    #                         # For a real analysis we would use V = 10.
    #                         V = 10,
    #                         SL.library = sl.lib)
    #
    # # We run summary on the cv_sl object rather than simply printing the object.
    # summary(cv_sl)
    # # Review the distribution of the best single learner as external CV folds.
    # table(simplify2array(cv_sl$whichDiscreteSL))
    # plot(cv_sl) + theme_bw()

    # dont do this for now
    # mod_sd <- sd(exposure - mean_predict)
    # feature_names <- GPS_mod$feature_names
    # GPS <- dnorm(exposure, mean = mean_predict, sd = mod_sd) # probability of exposure given confounders
    # Nm <- dnorm(exposure, mean = mean(exposure, na.rm = TRUE), sd = sd(exposure, na.rm = TRUE)) # General probability of exposure
    # data_example$IPW <- Nm / GPS #stabilized propensity score for really low or high weights
    # 
    # Subset model to middle 95 of weights (drop observations with extreme weights and estimate with remining data + weights)
    # upper_bound = quantile(data_example$IPW, 0.975)
    lower_bound = 0
    upper_bound = quantile(W1$weights, 0.99)
    
    # Chop those weights down to
    data_example[IPW > upper_bound, IPW := upper_bound]
    data_ipw <- data_example
    data_ipw <- data_example
    #data_ipw = filter(data_example, between(IPW, lower_bound, upper_bound))
    
    # Make correlation table.... might be difficult 
    correlation_table <- 
        rbindlist(lapply(names(data.frame(confounders)), function(cov){
          pre_cor <- abs(cor(data_ipw$exposure, data_ipw[[cov]]))
          post_cor <- abs(wtd.cor(data_ipw$exposure, data_ipw[[cov]], data_ipw$IPW))
          data.table(covariate = cov, pre_cor = pre_cor, post_cor = post_cor[1])
        }))
    
    # abs_cor <- 
    #   correlation_table %>% 
    #   pivot_longer(-covariate) %>% 
    #   group_by(name) %>% 
    #   dplyr::summarize(abs_cor = mean(abs(value)))
    # 
    # gg_correlation <- 
    #   correlation_table %>% 
    #   pivot_longer(-covariate) %>% 
    #   mutate(name = factor(name, levels = unique(name))) %>% 
    #   ggplot(aes(x = covariate, y = value, color = name, group = name)) + 
    #   geom_point() +
    #   geom_line() + 
    #   theme_bw() +
    #   coord_flip()
    
    # Dropping these models now 
    
    propensity_lm <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_ipw, family = family, weights = IPW)
    #propensity_lm_nocf <- glm(Y ~ exposure, data = data_ipw, family = family, weights = IPW)
    
    propensity_gam <- mgcv::gam(Y ~ s(exposure) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_ipw, family = family, weights = IPW)
    #propensity_gam_nocf <- mgcv::gam(Y ~ s(exposure), data = data_ipw, family = family, weights = IPW)
    
    # Now try fitting the kennedy double robust method (errors right now)
    # cts_dr <- function(a, y, x, a.vals = seq(min(a), max(a), length.out = 100),
    #                    span = NULL, span.seq = seq(0.15, 1, by = 0.05), k = 10, # select span for loess
    #                    sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.glmnet", "SL.ranger"))
    # 
    # kennedy <- cts_dr(a = data_example$exposure, y = data_example$Y, x = data_example[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")])
  }
  

  
  # Now many points to evaluate the integral
  discrete_points = 500
  
  # Add predictions as new column
  # data_prediction  <- data_example[, .(exposure = seq(min(data_example$exposure),
  #                                                     max(data_example$exposure),
  #                                                     length.out = discrete_points),
  #                                      cf1 = 0,
  #                                      cf2 = 0,
  #                                      cf3 = 0,
  #                                      cf4 = 0,
  #                                      cf5 = 0,
  #                                      cf6 = 0)]
  # 
  data_prediction  <- data_example[, .(exposure = seq(0, 20,length.out = discrete_points),
                                       cf1 = 0,
                                       cf2 = 0,
                                       cf3 = 0,
                                       cf4 = 0,
                                       cf5 = 0,
                                       cf6 = 0)]
  
  # Get true value for outcome
  if (relationship == "linear") {
    if (family == "gaussian") {
      y_true <- as.numeric(20 + 0.1 *  data_prediction$exposure)
    } else if (family == "poisson")
      y_true <- exp(2 + 0.1 *  data_prediction$exposure)
  } else if (relationship == "sublinear") {
    if (family == "gaussian") {
      # y_true <- as.numeric(20 - c(2, 2, 3, -1) %*% t(as.matrix(data_prediction[, .(cf1, cf2, cf3, cf4)])) - 2 * data_prediction$cf5 - 2 * data_prediction$cf6 + 5 * sqrt(data_prediction$exposure))
      y_true <- as.numeric(20 + 8 * log10(data_prediction$exposure + 1))
    } else if (family == "poisson") {
      y_true <- exp(2 + 1.5 * log10(data_prediction$exposure + 1))
    }
  } else if (relationship == "threshold") {
    if (family == "gaussian") {
      thresh_exp <- data_prediction$exposure
      thresh_exp[thresh_exp <= 5] <- 0
      thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
      y_true <- as.numeric(20 + 1 * thresh_exp)
    } else if (family == "poisson") {
      thresh_exp <- data_prediction$exposure
      thresh_exp[thresh_exp <= 5] <- 0
      thresh_exp[thresh_exp > 5] <- thresh_exp[thresh_exp > 5] - 5
      y_true <- exp(2 + 0.1 * thresh_exp)
    }
  } else {
    stop("Must either be linear, sublinear, or threshold as of now")
  }
  
  
  
  # Add predictions onto model
  if (adjust_confounder == T) {
    data_prediction[order(exposure), `:=`(linear_model = predict(linear_fit, data_prediction, type = "response"),
                                          gam_model = predict(gam_fit, newdata = data_prediction, type = "response"),
                                          linear_gps = predict(propensity_lm, newdata = data_prediction, type = "response"),
                                          gam_gps = predict(propensity_gam, newdata = data_prediction, type = "response"),
                                          true_fit = y_true)]
  } else {
    data_prediction[order(exposure), `:=`(linear_model = predict(linear_fit, data_prediction, type = "response"),
                                          gam_model = predict(gam_fit, newdata = data_prediction, type = "response"),
                                          true_fit = y_true)]
  }
  # Calculate trimmed reference exposure value
  trimmed_reference <- sort(data_prediction$exposure)[1]
  
  # Remove influence of the intercept to just compare risk difference or relative risk 
  if (family == "gaussian") {
    if (adjust_confounder) {
    model_types <- c("linear_model", "gam_model", "linear_gps", "gam_gps")
    data_prediction <- 
      data_prediction %>% 
      filter(exposure >= trimmed_reference) %>% 
      mutate(linear_model = linear_model - data_prediction[exposure == trimmed_reference, linear_model],
             gam_model = gam_model - data_prediction[exposure == trimmed_reference, gam_model],
             linear_gps = linear_gps - data_prediction[exposure == trimmed_reference, linear_gps],
             gam_gps = gam_gps - data_prediction[exposure == trimmed_reference, gam_gps],
             true_fit = true_fit - data_prediction[exposure == trimmed_reference, true_fit])
    } else {
      model_types <- c("linear_model", "gam_model")
      data_prediction <- 
        data_prediction %>% 
        filter(exposure >= trimmed_reference) %>% 
        mutate(linear_model = linear_model - data_prediction[exposure == trimmed_reference, linear_model],
               gam_model = gam_model - data_prediction[exposure == trimmed_reference, gam_model],
               true_fit = true_fit - data_prediction[exposure == trimmed_reference, true_fit])
    }
    
    # old way of doing it with using min as reference
    # data_prediction <- 
    #   data_prediction %>% 
    #   mutate(linear_model = linear_model - quantile(data_prediction$linear_model, 0.025),
    #          gam_model = gam_model - data_prediction$gam_model[1],
    #          true_fit = true_fit - data_prediction$true_fit[1])
    # 
    
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
  
  # Return plot if asked
  if (just_plots == T) {
    plot_fits <- 
      data_prediction %>% 
      pivot_longer(c(all_of(model_types), "true_fit"), names_to = "model", values_to = "prediction") %>%
      arrange(exposure) %>%
      ggplot(aes(x = exposure, y = prediction, color = model, linetype = model)) +
      geom_line() + 
      labs(x = "Exposure concentration", y = "Relative risk of death")
    
    return(plot_fits)
  }
  
  data_metrics <- 
    data_prediction %>% 
    pivot_longer(model_types, names_to = "model", values_to = "prediction") %>% 
    data.table()
  
  data_metrics <- 
    data_metrics[, .(bias = prediction - true_fit,
                     mse = (prediction - true_fit) ^ 2), by = .(model)]
  
  data_metrics <- data_metrics[, .(bias = mean(bias), mse = mean(mse)), by = .(model)]
  return_columns <- c("exposure", model_types, "true_fit")
  return(list(metrics = data_metrics, predictions = data_prediction[, ..return_columns], cor_table = correlation_table))
}


# Function for evaluating eSCHIF 
## Currently that will take forever to run because eSCHIF runs 1,000 versions of its funciton anyway, so thinking if that is worth it. 




