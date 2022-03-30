wrapper <- function(data, a.vals, n.boot = 1000,
                    sex = c("both", "male", "female"),
                    race = c("all", "white", "black", "hispanic", "asian")) {
  
  zip_cov <- c("pm25", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "poverty", "education",
               "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax", "region") 
  ind_cov <- c("dual", "entry_age_break", "followup_year")
  
  if (sex == "male") {
    sex <- 1
  } else if (sex == "female") {
    sex <- 2
  } else {
    sex <- c(1,2)
  }
  
  if (race == "white") {
    race <- 1
  } else if (race == "black") {
    race <- 2
  } else if (race == "hispanic") {
    race <- 5
  } else if (race == "asian") {
    race <- 4
  } else {
    race <- c(0,1,2,3,4,5,6)
  }
  
  sub_data <- data %>% filter(data$race %in% race & data$sex %in% sex)
  
  # Covariates and Outcomes
  outcomes <- data.table(zip = sub_data$zip, year = sub_data$year,
                         dead = sub_data$dead, time_count = sub_data$time_count)[,lapply(.SD, sum), by = c("zip", "year")]
  
  cmat <- data.frame(model.matrix(~ ., data = sub_data[,ind_cov])[,-1])
  ind_covariates <- data.frame(zip = outcomes$zip, year = outcomes$year)
  
  for (i in 1:ncol(cmat)){
    
    mat <- data.table(zip = sub_data$zip, year = sub_data$year, time_count = sub_data$time_count, val = cmat[,i])
    ind_tmp <- mat[,list(weighted.mean(val,time_count)), by = c("zip", "year")]
    colnames(ind_tmp)[3] <- colnames(cmat)[i]
    ind_covariates <- merge(ind_covariates, ind_tmp, by = c("zip", "year"), all.x = TRUE)
    
  }
  
  zip_covariates <- data.table(zip = sub_data$zip, year = sub_data$year,
                               model.matrix(~ ., data = sub_data[,zip_cov])[,-1])[,lapply(.SD, min), by = c("zip", "year")]
  
  new_data <- merge(outcomes, ind_covariates, by = c("zip", "year")) %>%
    merge(zip_covariates, by = c("zip", "year"))
  new_data$zip <- factor(new_data$zip)
  new_data$year <- factor(new_data$year)
  new_data.list <- split(data, list(data$zip))
  n_zip <- length(unique(data$zip))
  
  rm(data, sub_data, outcomes, zip_covariates, ind_covariates, 
     ind_tmp, cmat, mat, race, sex, ind_cov, zip_cov); gc()
  
  out <- mclapply(1:(n.boot + 1), mc.cores = 24, function(j, new_data, new_data.list, a.vals) {
    
    if (j == 1){
      boot_data <- new_data
    } else{  
      idx <- sample(1:n_zip, floor(2*sqrt(n_zip)), replace = TRUE) 
      boot_data <- data.frame(Reduce(rbind, new_data.list[idx]))
    }
    
    x <- subset(boot_data, select = -c(zip, dead, time_count, pm25))
    a <- boot_data$pm25
    y <- boot_data$dead
    offset <- log(boot_data$time_count)
    
    estimate <- fit_models(a = a, x = x, y = y, offset = offset, a.vals = a.vals)
    
    return(estimate)
    
  }, new_data = new_data, new_data.list = new_data.list, a.vals = a.vals)
  
  dr_estimate <- out[[1]]
  dr_boot <- out[2:(n.boot + 1)]
  out <- data.frame(a.vals = a.vals, estimate = dr_estimate, Reduce(rbind, dr_boot))
  colnames(out) <- c("a.vals", "estimate", paste0("boot", 1:n.boot))
  
  return(out)
  
}

fit_models <- function(a, x, y, offset, a.vals) {
  
  ## Matching Estimator
  # match_pop <- generate_pseudo_pop(Y = y, w = a, c = x,
  #                                  ci_appr = "matching", pred_model = "sl",
  #                                  gps_model = "parametric", use_cov_transfoqd = TRUE,
  #                                  transfoqders = list("pow2", "pow3"), sl_lib = c("m_xgboost"),
  #                                  params = list(xgb_nrounds=c(50)), nthread = 8, # number of cores, you can change,
  #                                  covar_bl_method = "absolute", covar_bl_trs = 0.1,
  #                                  trim_quantiles = c(0.025,0.975), # trimmed, you can change,
  #                                  optimized_compile = TRUE, #created a column counter for how many times matched,
  #                                  max_attempt = 5, matching_fun = "matching_l1",
  #                                  delta_n = (max(a.vals) - min(a.vals)), scale = 1.0)
  # 
  # match_data <- match_pop$pseudo_pop
  # match_data$offset <- offset[match_data$row_index]
  # match_data <- subset(match_data, counter > 0)
  # match_curve <- mgcv::bam(Y ~ s(w, bs = 'tp'), data = match_data, offset = match_data$offset, 
  #                          family = poisson(link = "log"), weights = match_data$counter)
  
  ## Doubly-Robust Estimator
  n <- length(a)
  x <- data.frame(x)
  
  # set up evaluation points & matrices for predictions
  xa.new <- rbind(cbind(x, a = a), cbind(x[rep(1:n, length(a.vals)), ], a = rep(a.vals, each = n)))
  xa.new <- data.frame(xa.new)
  colnames(x) <- colnames(xa.new)[-ncol(xa.new)]
  
  # estimate nuisance functions via super learner
  pimod <- ranger(a ~ ., data = xa.new[1:n,], num.trees = 200, min.node.size = 50)
  pimod.vals <- predict(pimod, data = xa.new[,-ncol(xa.new)], type = "response")$predictions
  pi2mod <- mean((a - pimod.vals[1:n])^2)
  fmla <- formula(paste("y ~ 0 + ns(a, df = 5, intercept = TRUE) +", paste(colnames(x), collapse = "+")))
  mumod <- glm(fmla, data = xa.new[1:n,], offset = offset, family = poisson(link = "log"))
  muhat.vals <- predict(mumod, newdata = xa.new, type = "response")
  
  # construct estimated pi/varpi and mu/m values
  a.std <- (xa.new$a - pimod.vals) / sqrt(pi2mod)
  pihat.vals <- approx(density(a.std[1:n])$x, 
                       density(a.std[1:n])$y, 
                       xout = a.std)$y / sqrt(pi2mod)
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  
  # truncation
  pihat[pihat < quantile(pihat, 0.025)] <- quantile(pihat, 0.025)
  pihat[pihat > quantile(pihat, 0.975)] <- quantile(pihat, 0.975)
  phat[phat < quantile(phat, 0.025)] <- quantile(phat, 0.025)
  phat[phat > quantile(phat, 0.975)] <- quantile(phat, 0.975)
  
  # tmle steps
  nsa <- ns(a.std[1:n], df = 5, intercept = TRUE)
  g <- nsa*phat/pihat
  mustard <- glm(y ~ 0 + ., data = data.frame(y = y, g), 
                 offset = log(muhat) + offset, family = poisson(link = "log"))
  
  phat_all <- rep(colMeans(pihat.mat), each = n)
  h <- predict(nsa, new.x = a.std[-(1:n)])*phat_all/pihat.vals[-(1:n)]
  epsilon <- coefs(mustard)
  est.vals <- exp(log(muhat.vals[(-1:n)]) + c(h%*%epsilon))
  est.mat <- matrix(est.vals, nrow = n, ncol = length(a.vals))
  estimate <- colMeans(est.mat)
  
  return(estimate)
  
}