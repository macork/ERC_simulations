# Xiao's old matching functions that I reference when understanding causalGPS
#' GPS Matching function to create matched set.
#'
#' @param Y a vector of observed outcome variable.
#' @param w a vector of observed continuous exposure variable.
#' @param c a data frame or matrix of observed covariates variable.
#' @param matching_fun a specifed matching function (Default is "matching_l1" (Manhattan distance matching)).
#' @param scale a specified scale parameter to control the relative weight that is attributed to the distance measures of the exposure versus the GPS estimates (Default is 0.5).
#' @param delta_n a specified caliper parameter on the exposure (Default is 1).
#' @param sl.lib a set of methods used for estimating GPS (Default is ("SL.xgboost","SL.earth","SL.gam","SL.ranger")).
#' @return
#' \code{matched_set}: The function returns a data.table saved the constructed matched set by the proposed GPS matching approaches.
#' @export

# Create matched set using GPS matching approaches
create_matching <- function(Y,
                            w,
                            c,
                            matching_fun = matching_l1,
                            sl.lib = c("SL.xgboost", "SL.earth", "SL.gam", "SL.ranger"),
                            scale = 0.5,
                            delta_n = 1){
  ## GPS function estimation
  e_gps <- SuperLearner(Y = w, X = data.frame(c), SL.library = sl.lib)
  e_gps_pred <- e_gps$SL.predict
  e_gps_std <- SuperLearner(Y = abs(w - e_gps_pred), X = c, SL.library = sl.lib)
  e_gps_std_pred <- e_gps_std$SL.predict
  w_resid <- (w-e_gps_pred) / e_gps_std_pred
  gps <- approx(density(w_resid, na.rm = TRUE)$x,
                density(w_resid, na.rm = TRUE)$y,
                xout = w_resid,
                rule = 2)$y
  
  dataset <- cbind(Y, w, c, gps)
  
  bin.num <- seq(min(dataset$w) + delta_n / 2, max(dataset$w), by = delta_n)
  
  matched_set <-  lapply(bin.num,
                         matching_fun,
                         dataset = dataset,
                         e_gps_pred = e_gps_pred,
                         e_gps_std_pred = e_gps_std_pred,
                         delta_n = delta_n,
                         w_resid = w_resid,
                         scale = scale)
  
  return(data.table(Reduce(rbind, matched_set)))
}

# Matching function using L1 distance on single exposure level w
matching_l1 <- function(dataset,
                        e_gps_pred,
                        e_gps_std_pred,
                        w,
                        delta_n = 1,
                        w_resid,
                        scale)
{
  w_new <- (w - e_gps_pred) / e_gps_std_pred
  p.w <- approx(density(w_resid, na.rm = TRUE)$x,
                density(w_resid, na.rm = TRUE)$y,
                xout = w_new,
                rule = 2)$y
  
  w.min <- min(dataset[["w"]], na.rm = T)
  w.max <- max(dataset[["w"]], na.rm = T)
  gps.min <- min(dataset[["gps"]], na.rm = T)
  gps.max <- max(dataset[["gps"]], na.rm = T)
  ##
  dataset <- transform(dataset,
                       std.w = (w - w.min) / (w.max - w.min),
                       std.gps = (gps - gps.min) / (gps.max - gps.min))
  std.w <- (w - w.min) / (w.max - w.min)
  std.p.w <- (p.w - gps.min) / (gps.max - gps.min)
  ##
  dataset.subset <- dataset[abs(dataset[["w"]] - w) <= (delta_n / 2), ]
  ##
  wm <- apply(abs(outer(dataset.subset[["std.gps"]], std.p.w, `-`)) * scale,
              2,
              function(x) which.min(abs(dataset.subset[["std.w"]] - std.w) * (1 - scale) + x)
  )
  dp <- dataset.subset[wm, ]
  return(dp)
  gc()
}

absolute_corr_fun <- function(w,
                              c){
  absolute_corr<- sapply(colnames(c),function(i){
    abs(cor(w,c[[i]],method = c("spearman")))})
  
  return(list(absolute_corr = absolute_corr,
              mean_absolute_corr = mean(absolute_corr)))
}


#' Estimate smoothed exposure-response function (ERF).
#'
#' @param matched_Y a vector of outcome variable in matched set.
#' @param matched_w a vector of continuous exposure variable in matched set.
#' @param bw.seq a vector of bandwidth values (Default is seq(0.2,2,0.2)).
#' @param w.vals a vector of exposure levels that ERF curves was evaluated.
#' @return
#' \code{erf}: The function returns a vector saved the output values of exposure-response function (ERF) given input \code{w.vals}.
#' @export

# Fit non-parametric kernel smoothing on matched set
matching_smooth<-function(matched_Y,
                          matched_w,
                          bw.seq = seq(0.2, 2, 0.2),
                          w.vals){
  ## The specified Gaussian kernel
  kern_fun <- function(t){ dnorm(t) }
  w_fun <- function(bw){
    w.avals <- NULL
    for (w.val in w.vals){
      w.std <- (matched_w-w.val) / bw
      kern.std <- kern_fun(w.std) / bw
      w.avals <- c(w.avals, mean(w.std^2*kern.std)*(kern_fun(0)/bw) /
                     (mean(kern.std) * mean(w.std^2 * kern.std) - mean(w.std * kern.std)^2))
    }
    return(w.avals / length(matched_w))
  }
  hatvals <- function(bw){approx(w.vals,
                                 w_fun(bw),
                                 xout = matched_w,
                                 rule = 2)$y}
  smooth_fun <- function(out,bw){
    approx(locpoly(matched_w,
                   out,bandwidth = bw,
                   gridsize = 1000),
           xout = matched_w, rule = 2)$y
  }
  ##
  risk_fun <- function(h){
    hats <- hatvals(h); mean(((matched_Y - smooth_fun(matched_Y, bw = h)) / (1 - hats))^2)
  }
  risk.val <- sapply(bw.seq, risk_fun)
  h.opt <- bw.seq[which.min(risk.val)]
  
  erf <- approx(locpoly(matched_w, matched_Y, bandwidth = h.opt), xout = w.vals)$y
  return(erf)
}


F.aac.iter=function(i,data,ps.model,ps.num,rep,criterion) {
  # i: number of iterations (trees)
  # data: dataset containing the treatment and the covariates
  # ps.model: the boosting model to estimate p(T_i|X_i)
  # ps.num: the estimated p(T_i)
  # rep: number of replications in bootstrap
  # criterion: the correlation metric used as the stopping criterion
  GBM.fitted=predict(ps.model,newdata=data,n.trees=floor(i),
                     type="response")
  ps.den=dnorm((data$T-GBM.fitted)/sd(data$T-GBM.fitted),0,1)
  wt=ps.num/ps.den
  aac_iter=rep(NA,rep)
  for (i in 1:rep){
    bo=sample(1:dim(data)[1],replace=TRUE,prob=wt)
    newsample=data[bo,]
    j.drop=match(c("T"),names(data))
    j.drop=j.drop[!is.na(j.drop)]
    x=newsample[,-j.drop]
    if(criterion=="spearman"|criterion=="kendall"){
      ac=apply(x, MARGIN=2, FUN=cor, y=newsample$T,
               method=criterion)
    } else if (criterion=="distance"){
      ac=apply(x, MARGIN=2, FUN=dcor, y=newsample$T)
    } else if (criterion=="pearson"){
      ac=matrix(NA,dim(x)[2],1)
      for (j in 1:dim(x)[2]){
        ac[j]=ifelse (!is.factor(x[,j]), cor(newsample$T, x[,j],
                                             method=criterion),polyserial(newsample$T, x[,j]))
      }
    } else print("The criterion is not correctly specified")
    aac_iter[i]=mean(abs(1/2*log((1+ac)/(1-ac))),na.rm=TRUE)
  }
  aac=mean(aac_iter)
  return(aac)
}

x = data_gps[, c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]
w = data_gps$exposure

library(gbm)
library(polycor)
mydata = data.frame(T=w,X=x)
model.num = lm(T~1,data=mydata)
ps.num=dnorm((mydata$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
model.den=gbm(T~.,data=mydata, shrinkage = 0.0005,
              interaction.depth=4, distribution="gaussian",n.trees=20000)
opt=optimize(F.aac.iter,interval=c(1,20000), data=mydata, ps.model=model.den,
             ps.num=ps.num,rep=50,criterion="pearson")
best.aac.iter=opt$minimum
best.aac=opt$objective


model.den$fitted=predict(model.den,newdata=mydata,n.trees=floor(best.aac.iter), type="response")
ps.den=dnorm((mydata$T-model.den$fitted)/sd(mydata$T-model.den$fitted),0,1)
weight.gbm=ps.num/ps.den

data_gps$IPW <- weight.gbm


