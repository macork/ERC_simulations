
# wrapper function to fit a doubly-robust ERC using LOESS regression on nonparametric models
cts_dr <- function(a, y, x, a.vals = seq(min(a), max(a), length.out = 100), 
                   span = NULL, span.seq = seq(0.15, 1, by = 0.05), k = 10, # select span for loess
                   sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.glmnet", "SL.ranger")) {	
  
  n <- length(y)
  
  wrap <- np_est(y = y, a = a, x = x, sl.lib = sl.lib)
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  pihat <- wrap$pihat
  phat <- wrap$phat
  int <- wrap$int
  wts <- wrap$data$wts
  psi <- (y - muhat)/(pihat/phat) + mhat
  
  if(is.null(span)) {
    # What size should M be? is this n?
    folds <- sample(x = k, size = n, replace = TRUE)
    
    cv.mat <- sapply(span.seq, function(h, ...) {
      
      cv.vec <- rep(NA, k)
      
      for(j in 1:k) {
        
        # Running into errors here?
        preds <- sapply(a[folds == j], dr_est, psi = psi[folds != j], a = a[folds != j], 
                        int = int[folds != j], wts = wts[folds != j], span = h, se.fit = FALSE)
        cv.vec[j] <- mean(wts[folds == j]*(psi[folds == j] - preds)^2, na.rm = TRUE)
        # some predictions result in `NA` because of the `x` ranges in each fold
        
      }
      
      return(cv.vec)
      
    })
    
    cv.err <- colMeans(cv.mat)
    span <- span.seq[which.min(cv.err)]
    
  }
  
  dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = int, wts = wts, span = span, se.fit = TRUE)
  
  estimate <- dr_out[1,]
  variance <- dr_out[2,]
  
  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance)	
  
  return(out)
  
}

# LOESS function
dr_est <- function(newa, a, psi, int, wts, span, se.fit = FALSE) {
  
  a.std <- a - newa
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  knn <- rep(0, length(a))
  knn[idx] <- 1
  max.a.std <- max(abs(a.std*knn))
  k.std <- wts*((70/81)*(1 - abs(a.std/max.a.std)^3)^3)*knn #tricube distance
  
  ## Gaussian kernel smoothing
  # a.std <- (a - newa) / h
  # k.std <- dnorm(a.std) / h
  
  gh <- cbind(1, a.std)
  gh.inv <- solve(t(gh) %*% diag(k.std) %*% gh)
  b <- optim(par = c(0,0), fn = opt_fun, k.std = k.std, psi = psi, gh = gh)
  mu <- b$par[1]
  
  if (se.fit){
    
    v.inf <- (psi + int - c(gh%*%b$par))^2
    sig <- gh.inv %*% t(gh) %*% diag(k.std) %*% diag(v.inf) %*% diag(k.std) %*% gh %*% gh.inv
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)
  
}

# Nonparametric nuisance parameter estimation
np_est <- function(a, y, x, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.glmnet", "SL.ranger")){
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a)
  
  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = y, X = xa, SL.library = sl.lib, family = gaussian())
  mumod.vals <- c(mumod$SL.predict)
  muhat <- c(mumod$SL.predict)
  # Add this here
  #muhat = predict(mumod) or mumod$fitted.values
  
  
  # estimate nuisance GPS functions via super learner
  pimod <- SuperLearner(Y = a, X = x, SL.library = sl.lib, family = gaussian())
  pimod.vals <- c(pimod$SL.predict)
  # model variance
  pi2mod <- SuperLearner(Y = c(a - pimod.vals)^2, X = x, SL.library = "SL.ranger", family = gaussian()) # you can make this faster/fancier
  pi2mod.vals <- c(pi2mod$SL.predict)
  
  # exposure models
  pihat <- dnorm(a, pimod.vals, sqrt(pi2mod.vals))
  phat <- sapply(a, function(a.tmp, ...) mean(dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals))))
  
  # Makes assumption outcome is normal 
  # predict marginal outcomes given a.vals (or a)
  # muhat.mat <- sapply(a2, function(a.tmp, ...) {
  #   
  #   xa.tmp <- data.frame(x, a = a.tmp)
  #   colnames(xa.tmp) <- colnames(xa) 
  #   return(c(predict(mumod, newdata = xa.tmp)$pred))
  #   
  # })
  muhat.mat <- 
    do.call(cbind, mclapply(a, mc.cores = 10, function(a.tmp, ...) {
      
      xa.tmp <- data.frame(x, a = a.tmp)
      colnames(xa.tmp) <- colnames(xa) 
      return(c(predict(mumod, newdata = xa.tmp)$pred))
      
    }))
  # empirical integral from kennedy paper influence function
  # marginalize muhat.mat and integrate for influence curve
  mhat <- colMeans(muhat.mat)
  mhat.mat <- matrix(rep(mhat, n), byrow = TRUE, nrow = n)
  int <- rowMeans(muhat.mat - mhat.mat)
  
  out <- list(muhat = muhat, mhat = mhat, pihat = pihat, phat = phat, int = int)
  
  return(out)
  
}

# optimization used in dr_est
opt_fun <- function(par, k.std, psi, gh) {
  
  sum(k.std*(psi - c(gh %*% par))^2)
  
}