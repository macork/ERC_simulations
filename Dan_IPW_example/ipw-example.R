# This code was created by Dan as an example to show where IPW weighting was improving effect estimates as compared to fitting a linear model
# Most of this has been moved to simulations_confounding function 


rm(list = ls())
library(tidyverse)
library(data.table)
library(MASS)
library(ggeffects)
library(parallel)
library(xgboost)
library(SuperLearner)
library(weights)
library(CBPS)
library(caret)
library(WeightIt)
library(chngpt)
library(cobalt)
library(KernSmooth)
library("devtools")
#install_github("fasrc/CausalGPS")
library("CausalGPS")


n <- 1000
res <- list(W = list(), DR = list())
cor_ipw <- list()
cor_ipw2 <- list()
cor_ent <- list()


set.seed(23)
cf <- mvrnorm(n = sample_size,
              mu = rep(0, 4),
              Sigma = diag(4))
cf5 <- sample(c((-2):2), sample_size, replace = T)
cf6 <- runif(sample_size, min = -3, max = 3)
cov = data.frame(cbind(cf, cf5, cf6))
colnames(cov) = c("c1", "c2", "c3", "c4", "c5", "c6")

scale_exposure <- function(x) {
  20 * (x - min(x)) / (max(x) - min(x))
}


w = 10 * cov_function(cov) + rnorm(sample_size, mean = 0, sd = 1)
w = scale_exposure(w)
y = as.numeric(20 - c(2, 2, 3, -1, -2, -2) %*% t(cov) + 0.1 * w + rnorm(sample_size, mean = 0, sd = 10))

# Simulate covariates, exposure (w) and outcome (y) as interaction of all three
#cov <- data.frame(c1 = runif(n), c2 = rnorm(n), c3 = rbinom(n, 1, 0.5))
#w <- rnorm(n, cov$c1 + cov$c2 * cov$c3)
#y <- 0.1 * w + cov$c1 + cov$c2 * cov$c3 + rnorm(n)

# Train xgboost model
tune_control <- 
  caret::trainControl(method = "cv",
                      number = 5,
                      verboseIter = FALSE,
                      allowParallel = TRUE)
nrounds <- 400
tune_grid <- expand.grid(nrounds = seq(from = 50, to = nrounds, by = 50),
                         eta = c(0.025, 0.05, 0.1, 0.15),
                         max_depth = c(2),
                         gamma = 0,
                         colsample_bytree = 1,
                         min_child_weight = c(1,2),
                         subsample = 1)
xgbt_model <- 
  train(w ~ ., data = data.frame(cbind(w, cov)), method = "xgbTree",
        trControl = tune_control, tuneGrid = tune_grid, verbosity = 0)

# Get predictions from model and compare
predVals <- extractPrediction(list(xgbt_model))
plotObsVsPred(predVals)

# Create weighted IPW estimator in denominator 
ipw2 <- dnorm(w, mean(w), sd(w)) / dnorm(w, predVals$pred, sqrt(mean(resid(xgbt_model)^2)))
#ipw <- 1 / dnorm(w, predVals$pred, sqrt(mean(resid(xgbt_model)^2)))



correlation_table <- 
  rbindlist(lapply(names(data.frame(cov)), function(c){
    pre_cor <- abs(cor(w, cov[[c]]))
    post_cor <- abs(wtd.cor(w, cov[[c]], ipw2))
    data.table(covariate = c, pre_cor = pre_cor, post_cor = post_cor[1])
  }))

# Now extract best tune values for gps prediction function
pseudo_pop <- generate_pseudo_pop(y,
                                  w,
                                  cov,
                                  ci_appr = "matching",
                                  pred_model = "sl",
                                  gps_model = "non-parametric",
                                  use_cov_transform = TRUE,
                                  transformers = list("pow2", "pow3", "abs", "scale"),
                                  trim_quantiles = c(0.01,0.99),
                                  optimized_compile = TRUE,
                                  sl_lib = c("m_xgboost"),
                                  params = list(xgb_max_depth = c(xgbt_model$bestTune$max_depth, xgbt_model$bestTune$max_depth + 1),
                                                xgb_rounds = c(xgbt_model$bestTune$nrounds - 100, xgbt_model$bestTune$nrounds, xgbt_model$bestTune$nrounds + 100),
                                                xgb_eta = c(xgbt_model$bestTune$eta, xgbt_model$bestTune$eta / 2)),
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  covar_bl_trs_type = "mean",
                                  max_attempt = 10,
                                  matching_fun = "matching_l1",
                                  delta_n = 1,
                                  scale = 1,
                                  nthread = 1)

# Semi parametric erf
gam_model <- gam::gam(formula = Y ~ w, 
                    family = "gaussian", data = data.frame(pseudo_pop$pseudo_pop), weights = counter)
plot(gam_model)


CausalGPS::generate_pseudo_pop()

estimate_gps_out <- data_with_gps
pseudo_pop <- compile_pseudo_pop(dataset = estimate_gps_out, 
                                 ci_appr = ci_appr, gps_model, bin_seq, nthread = nthread, 
                                 trim_quantiles = trim_quantiles, optimized_compile = optimized_compile, 
                                 ...)


# Compare the correlations here, notice there is no covariate balance achieved
correlation_table_causal <- 
  rbindlist(lapply(c("c1", "c2", "c3"), function(c){
    pre_cor <- abs(cor(data_with_gps$w, data_with_gps[, c]))
    post_cor <- abs(wtd.cor(data_with_gps$w, data_with_gps[, c]))
    data.table(covariate = c, pre_cor = pre_cor, post_cor = post_cor[1])
  }))


# Run through Xiao's functions
source("~/Desktop/Francesca_research/Simulation_studies/Xiao_functions.R")
Y = y
c = cov

# Create matched set
matched_set = create_matching(Y,
                              w,
                              c,
                              matching_fun = matching_l1,
                              sl.lib = c("SL.xgboost","SL.earth","SL.gam","SL.ranger"),
                              scale = 1,
                              delta_n=1)


# Now 

erf = matching_smooth(matched_Y = matched_set$Y,
                      matched_w = matched_set$w,
                      bw.seq = seq(0.2,2,0.2),
                      w.vals = seq(min(w),max(w),length.out = 100))



w = y
c = cov
e_gps <- train_it(target = w, input = c, pred_model, 
                  sl_lib_internal = sl_lib_internal)
e_gps_pred <- e_gps$SL.predict
e_gps_std_pred <- stats::sd(w - e_gps_pred)
w_resid <- compute_resid(w, e_gps_pred, e_gps_std_pred)
gps <- stats::dnorm(w, mean = e_gps_pred, sd = e_gps_std_pred)

data_with_gps <- cbind(y, w, cov, gps)

# Need to load function here 
pr_mdl <- SuperLearner::SuperLearner(Y=y, X=data.frame(cov),
                                     SL.library = sl_lib_internal)
return(pr_mdl)





# Run a couple of times to test---------------------------------------------------------
# Here is where we can run our results now
set.seed(3818)
for (s in 1:10) {
  cov <- data.frame(c1 = runif(n), c2 = rnorm(n), c3 = rbinom(n, 1, 0.5))
  w <- rnorm(n, cov$c1 + cov$c2 * cov$c3)
  y <- 0.1 * w + cov$c1 + cov$c2 * cov$c3 + rnorm(n)
  
  # Propensity score model (correctly specified)
  ps_mod <- lm(w ~ c1 + c2 * c3, data = cov)
  ipw <- dnorm(w, mean(w), sd(w)) / dnorm(w, ps_mod$fitted.values, sqrt(mean(ps_mod$residuals^2)))
  
  # Using weights as described by Dan 
  ent <- weightit(w ~ c1 + c2 + c3, data = cbind(cov, w), method = "ebal", stabilize = T)
  
  # Using other weights
  #np_balance <- weightit(w ~ c1 + c2 + c3, data = cbind(cov, w), method = "npcbps", stabilize = T)
  # Use nonlinear term
  # Train xgboost model
  tune_control <- 
    caret::trainControl(method = "cv",
                        number = 5,
                        verboseIter = FALSE,
                        allowParallel = TRUE)
  nrounds <- 400
  tune_grid <- expand.grid(nrounds = seq(from = 50, to = nrounds, by = 50),
                           eta = c(0.025, 0.05, 0.1),
                           max_depth = c(2),
                           gamma = 0,
                           colsample_bytree = 1,
                           min_child_weight = c(1,2),
                           subsample = 1)
  xgbt_model <- 
    train(w ~ ., data = data.frame(cbind(w, cov)), method = "xgbTree",
          trControl = tune_control, tuneGrid = tune_grid, verbose = 0)
  
  # Get predictions from model and compare
  predVals <- extractPrediction(list(xgbt_model))
  plotObsVsPred(predVals)
  
  # Now get IPW weights from xgboost
  ipw2 <- dnorm(w, mean(w), sd(w)) / dnorm(w, predVals$pred, sqrt(mean(resid(xgbt_model)^2)))
  
  correlation_table <- 
    rbindlist(lapply(names(data.frame(cov)), function(c){
      pre_cor <- abs(cor(w, cov[[c]]))
      post_cor <- abs(wtd.cor(w, cov[[c]], ipw))
      data.table(covariate = c, pre_cor = pre_cor, post_cor = post_cor[1])
    }))
  
  
  # Using CausalGPS package ----
  # Start using superlearner (realize that the absolute correlation is different than expected)
  # Grid search for values, what happens over this grid? How many values do you look over?
  # 0.1 through 2 how to search through grids, run through causal GPS 100 times. Then used that for the entire simulation
  pseudo_pop <- generate_pseudo_pop(y,
                                    w,
                                    cov,
                                    ci_appr = "matching",
                                    pred_model = "sl",
                                    gps_model = "non-parametric",
                                    transformers = list("pow2", "pow3", "abs", "scale"),
                                    use_cov_transform = TRUE,
                                    trim_quantiles = c(0.01,0.99),
                                    optimized_compile = TRUE,
                                    sl_lib = c("m_xgboost"),
                                    params = list(xgb_max_depth = c(xgbt_model$bestTune$max_depth, xgbt_model$bestTune$max_depth + 1),
                                                  xgb_rounds = c(xgbt_model$bestTune$nrounds - 100, xgbt_model$bestTune$nrounds, xgbt_model$bestTune$nrounds + 100),
                                                  xgb_eta = c(xgbt_model$bestTune$eta, xgbt_model$bestTune$eta / 2)),
                                    covar_bl_method = "absolute",
                                    covar_bl_trs = 0.1,
                                    covar_bl_trs_type = "mean",
                                    max_attempt = 10,
                                    matching_fun = "matching_l1",
                                    delta_n = 0.5,
                                    scale = 1,
                                    nthread = 1)
  
  
  # erf_obj <- estimate_npmetric_erf(matched_Y = pseudo_pop$pseudo_pop$Y,
  #                                  matched_w = pseudo_pop$pseudo_pop$w,
  #                                  matched_counter = pseudo_pop$pseudo_pop$counter,
  #                                  bw_seq = seq(0.1, 2, length.out = 10),
  #                                  w_vals = seq(min(w),max(w),length.out = 10),
  #                                  nthread = 1)
  # 
  # plot(erf_obj)
  # 
  # # Do semi parametric erf instead of parameteric here 
  # gam_model <- gam::gam(formula = Y ~ s(w), 
  #                       family = "gaussian", data = data.frame(pseudo_pop$pseudo_pop), weights = counter)
  # 
  # gam_model <- gam::gam(formula = Y ~ w, 
  #                       family = "gaussian", data = data.frame(pseudo_pop$pseudo_pop), weights = counter)
  # plot(gam_model)
  # 
  # # Plotting gam model is helpful here 
  # 
  # 
  # erf_obj <- estimate_semipmetric_erf(formula = Y ~ w, 
  #                                 family = "gaussian", 
  #                                 data = data.frame(pseudo_pop$pseudo_pop),
  #                                 ci_appr = "matching")
  # plot(erf_obj)
  # 
  
  cor_ipw[[s]] <- cov.wt(cbind.data.frame(w, cov), ipw, TRUE)$cor[,1]
  cor_ipw2[[s]] <- cov.wt(cbind.data.frame(w, cov), ipw2, TRUE)$cor[,1]
  cor_ent[[s]] <- cov.wt(cbind.data.frame(w, cov), ent$weights, TRUE)$cor[,1]
  
  # Model with W only
  res[["W"]][[s]] <- c(reg = lm(y ~ w)$coef[2],
                       ipw = lm(y ~ w, weights = ipw)$coef[2],
                       ipw2 = lm(y ~ w, weights = ipw2)$coef[2],
                       entropy = lm(y ~ w, weights = ent$weights)$coef[2],
                       causal_GPS = coef(gam::gam(formula = Y ~ w, family = "gaussian", data = data.frame(pseudo_pop$pseudo_pop), weights = counter))[2])
  
  # Double robust model (misspecified outcome)
  res[["DR"]][[s]] <- c(reg = lm(y ~ w + c1 + c2 + c3, data = cov)$coef[2],
                        ipw = lm(y ~ w + c1 + c2 + c3, data = cov, weights = ipw)$coef[2],
                        ipw2 = lm(y ~ w + c1 + c2 + c3, data = cov, weights = ipw2)$coef[2],
                        entropy = lm(y ~ w + c1 + c2 + c3, data = cov, weights = ent$weights)$coef[2],
                        causal_GPS = coef(gam::gam(formula = Y ~ w + c1 + c2 + c3, family = "gaussian", data = data.frame(pseudo_pop$pseudo_pop), weights = counter))[2])
}

# Create better boxplot
boxplot(do.call(rbind, res[["W"]]), main = "IPW weight", names = c("unweighted", "linear_ipw", "nonlinear_ipw", "balancing_ipw", "Causal_GPS")); abline(h = 0.1, col = "red")

boxplot(do.call(rbind, res[["DR"]]), main = "Double robust",  names = c("unweighted", "linear_ipw", "nonlinear_ipw", "balancing_ipw", "Causal_GPS")); abline(h = 0.1, col = "red")
boxplot(abs(do.call(rbind, cor_ipw)[, -1]), main = "Balance", ylab = "Abs corr"); abline(h = 0.1, col = "red")
boxplot(abs(do.call(rbind, cor_ipw2)[, -1]), main = "Balance", ylab = "Abs corr"); abline(h = 0.1, col = "red")
boxplot(abs(do.call(rbind, cor_ent)[, -1]), main = "Balance", ylab = "Abs corr"); abline(h = 0.1, col = "red")


set.seed(3818)
for (s in 1:100) {
  cov <- data.frame(c1 = runif(n), c2 = rnorm(n), c3 = rbinom(n, 1, 0.5))
  w <- rnorm(n, cov$c1 * cov$c2 * cov$c3)
  y <- 0.1 * w + cov$c1^2 + cov$c2 * cov$c3 + rnorm(n)
  
  # Propensity score model (correctly specified)
  ps_mod <- lm(w ~ c1 * c2 * c3, data = cov)
  ipw <- dnorm(w, mean(w), sd(w)) / dnorm(w, ps_mod$fitted.values, sqrt(mean(ps_mod$residuals^2)))
  
  # Using weights as described by Dan 
  ent <- weightit(w ~ c1 + c2 + c3, data = cbind(cov, w), method = "ebal", stabilize = T)
  
  
  
  cor[[s]] <- cov.wt(cbind.data.frame(w, cov), ent$weights, TRUE)$cor[,1]
  
  # Model with W only
  res[["W"]][[s]] <- c(reg = lm(y ~ w)$coef[2],
                       ipw = lm(y ~ w, weights = ipw)$coef[2],
                       entropy = lm(y ~ w, weights = ent$weights)$coef[2])
  
  # Double robust model (misspecified outcome)
  res[["DR"]][[s]] <- c(reg = lm(y ~ w + c1 + c2 + c3, data = cov)$coef[2],
                        ipw = lm(y ~ w + c1 + c2 + c3, data = cov, weights = ipw)$coef[2],
                        entropy = lm(y ~ w + c1 + c2 + c3, data = cov, weights = ent$weights)$coef[2])
}


boxplot(do.call(rbind, res[["W"]]), main = "W only"); abline(h = 0.1, col = "red")
boxplot(do.call(rbind, res[["DR"]]), main = "Double robust"); abline(h = 0.1, col = "red")
boxplot(abs(do.call(rbind, cor)[, -1]), main = "Balance", ylab = "Abs corr"); abline(h = 0.1, col = "red")