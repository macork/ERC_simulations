### Script for causal GPS integration, eventually to be moved to simulations r markdown file

rm(list = ls())
library(tidyverse)
library(data.table)
library(MASS)
library(ggeffects)
library(parallel)


set.seed(23)
source("~/Desktop/Francesca_research/Simulation_studies/simulation_functions.R")
sample_size = 1000
# Define exposure first
exposure <- sort(rnorm(sample_size, mean = 10, sd = 3))
exposure[exposure < 0] <- 0 # Truncate to be less than zero

# Delete in a bit
cf <- mvrnorm(n = sample_size,
              mu = rep(0, 4),
              Sigma = diag(4))
cf5 <- sample(c((-2):2), sample_size, replace = T)
cf6 <- runif(sample_size, min = -3, max = 3)
confounders = cbind(cf, cf5, cf6)
colnames(confounders) = c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")


sim_results <-
  mclapply(1:50, mc.cores = 10, function(i) {
    metrics_from_data(exposure = exposure, confounders = confounders, relationship = "linear", sample_size = sample_size, family = "poisson")
  })



#### First try data that was generated from causalGPS package itself
## Forced me to reload a lot of things to get this package to even work, but did finally load the package.
# Working through vignette 
library("devtools")
install_github("fasrc/CausalGPS")
library("CausalGPS")


# Do an R markdown of the synthetic dataset and do the code, see if the code isn't working
# If there is something wrong with package, they can correct it 

set.seed(422)
n <- 10000
mydata <- generate_syn_data(sample_size=n)
year <- sample(x=c("2001","2002","2003","2004","2005"),size = n, replace = TRUE)
region <- sample(x=c("North", "South", "East", "West"),size = n, replace = TRUE)
mydata$year <- as.factor(year)
mydata$region <- as.factor(region)
mydata$cf5 <- as.factor(mydata$cf5)

pseudo_pop <- generate_pseudo_pop(mydata$Y,
                                  mydata$treat,
                                  mydata[c("cf1","cf2","cf3","cf4","cf5","cf6","year","region")],
                                  ci_appr = "matching",
                                  pred_model = "sl",
                                  gps_model = "non-parametric",
                                  use_cov_transform = TRUE,
                                  transformers = list("pow2", "pow3", "abs", "scale"),
                                  trim_quantiles = c(0.01,0.99),
                                  optimized_compile = TRUE,
                                  sl_lib = c("m_xgboost"),
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  max_attempt = 4,
                                  matching_fun = "matching_l1",
                                  delta_n = 1,
                                  scale = 0.5,
                                  nthread = 1)

plot(psuedo_pop)

pseudo_pop2 <- generate_pseudo_pop(mydata$Y,
                                  mydata$treat,
                                  mydata[c("cf1","cf2","cf3","cf4","cf5","cf6","year","region")],
                                  ci_appr = "matching",
                                  pred_model = "sl",
                                  gps_model = "non-parametric",
                                  use_cov_transform = TRUE,
                                  transformers = list("pow2", "pow3", "abs", "scale"),
                                  trim_quantiles = c(0.01,0.99),
                                  optimized_compile = TRUE,
                                  sl_lib = c("m_xgboost"),
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  max_attempt = 4,
                                  matching_fun = "matching_l1",
                                  delta_n = 1.5,
                                  scale = 1,
                                  nthread = 1)   

plot(pseudo_pop)

# This is code for generating simulation data 
erf_obj <- estimate_npmetric_erf(pseudo_pop$pseudo_pop$Y,
                                 pseudo_pop$pseudo_pop$w,
                                 bw_seq=seq(0.2,2,0.2),
                                 w_vals = seq(2,20,0.5),
                                 nthread = 1)

plot(erf_obj)


# Generat GPS value (internal function?)
data_with_gps <- estimate_gps(Y = mydata$Y,
                              w = mydata$treat,
                              c = mydata[c("cf1","cf2","cf3","cf4","cf5","cf6","year","region")],
                              pred_model = "sl",
                              internal_use = FALSE,
                              params = list(xgb_max_depth = c(3,4,5),
                                            xgb_rounds = c(10,20,30,40)),
                              nthread = 1,                                
                              sl_lib = c("m_xgboost"))



# Assessing covariate balance is important here, that is the step where we make sure the covariates are properly balanced...
# Working with Jenny and Priyanka on the causalGPS format 
# Matching portion is to make sure that the distribution of pre exposure covariates is similar across both groups 

# Why do we do the matching this way? Is it a problem if many of the same observed outcomes are used? 
# We want aboslute correlation to fall below 0.1, becuase theoretically they should be about zero. Not sure where the 0.1 comes from, but its interesting 
# Use global measures of absolute correlation to find the hyperparameters that give the best covariate balance here 
# What happens if we don't get good correlation matching? 

# Different simulation settings, try those out as well 
Here is an interesting tidbit:
  We also show an example where our GPS matching approach is not able to achieve covariate balance under a scenario where the covariates are strongly associated with the exposure (see Section S.3.5 of the Supplementary Materials). We suggest that researchers proceed to the analysis stage only if covariate balance has been achieved in the design stage.




# Uses absolute correlation


# Could use base case of linear confounding for first simulation, discuss with Dan if that is a good use of my time? I wanted some knowledge of this algorithm before 
# Talking with PRiiyanka and Jenny about using this package 

W =9×γ(C)+17+N(0,5);


# We specify the location constants (17, 22, 15, −6, −18, 13) and scale constants (9, 15, 9, 49, 42, 7) to 
# ensure all simulation scenarios 1-6 generate an exposure W with [5%, 95%] quantiles at approximately [0, 20], and thus all simulation scenarios are comparable, and align with the exposure range in our data application.
