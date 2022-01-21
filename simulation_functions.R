# File that contains all the functions I am using for the simulation 

data_generate_a <-
  function(sample_size = 1000, exposure = NA, confounders = NA, relationship = NA, family = NA, outcome_sd = 10) {
    
    # Abbreviate for clarity 
    cf <- confounders
    
    if (relationship == "linear") {
      if (family == "gaussian") {
        Y = as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(cf) + 0.1 * exposure + rnorm(sample_size, mean = 0, sd = outcome_sd))
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
metrics_from_data <- function(just_plots = F, exposure = NA, relationship = "linear", sample_size = 1000, 
                              confounders = NA, family = "gaussian", eschif_draws = NULL) {
  
  # Fit data generating mechanism
  data_example <- data_generate_a(sample_size = sample_size, exposure = exposure, confounders = confounders, relationship = relationship, family = family)
  
  # Now fit a linear model
  linear_fit <- glm(Y ~ exposure + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
  
  # Now fit a GAM model 
  gam_fit <- mgcv::gam(Y ~ s(exposure) + cf1 + cf2 + cf3 + cf4 + cf5 + cf6, data = data_example, family = family)
  
  # Now many points to evaluate the integral
  discrete_points = 1000
  
  # Add predictions as new column
  data_prediction  <- data_example[, .(exposure = seq(min(data_example$exposure),
                                                      max(data_example$exposure),
                                                      length.out = discrete_points),
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
  data_prediction[order(exposure), `:=`(linear_model = predict(linear_fit, data_prediction, type = "response"),
                                        gam_model = predict(gam_fit, newdata = data_prediction, type = "response"),
                                        true_fit = y_true)]
  # Calculate trimmed reference exposure value
  #trimmed_reference <- sort(data_prediction$exposure)[5]
  
  # Remove influence of the intercept to just compare risk difference or relative risk 
  if (family == "gaussian") {
    model_types <- c("linear_model", "gam_model")
    data_prediction <- 
      data_prediction %>% 
      filter(exposure >= trimmed_reference) %>% 
      mutate(linear_model = linear_model - data_prediction[exposure == trimmed_reference, linear_model],
             gam_model = gam_model - data_prediction[exposure == trimmed_reference, gam_model],
             true_fit = true_fit - data_prediction[exposure == trimmed_reference, true_fit])
    
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
  
  data_metrics[, `:=`(bias = prediction - true_fit,
                      mse = (prediction - true_fit)^2), by = .(model)]
  
  data_metrics <- data_metrics[, .(bias = mean(bias), mse = mean(mse)), by = .(model)]
  return_columns <- c("exposure", model_types, "true_fit")
  return(list(metrics = data_metrics, predictions = data_prediction[, ..return_columns]))
}


# Function for evaluating eSCHIF 
## Currently that will take forever to run because eSCHIF runs 1,000 versions of its funciton anyway, so thinking if that is worth it. 




