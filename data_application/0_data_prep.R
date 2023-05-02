# Code to prep data for model fitting (only need to run once for each new decision)
rm(list = ls())
library(data.table)
library(tidyverse)
library(chngpt)
library(WeightIt)
library(fastDummies)
library(mgcv)
library(Rcpp)
library(RcppEigen)
library("CausalGPS", lib.loc = "/n/home_fasse/mcork/apps/ERC_simulation/R_4.0.5")

# Name for input data
input_flag <- "kevin_trim_90"

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
out_dir <- paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/")
dir.create(out_dir)

# Load data from Kevin's workspace
load("/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/aggregate_data.RData")
data <- data.table(aggregate_data)

# See percent of age in each category
count(data, entry_age_break) %>% mutate(n = n / sum(n))

# See percent of each race in category
count(data, race) %>% mutate(n = n / sum(n))

# Calculate the mortality rate
data[, mort_rate := dead / time_count]

# now trim upper and lower quantiles for data analysis
data_trim <- data %>% select(all_of(c("zip", "pm25", "year"))) %>% unique()
lower_pm <- quantile(data_trim$pm25, 0.05)
upper_pm <- quantile(data_trim$pm25, 0.95)

data <-
  data %>%
  filter(between(pm25, lower_pm, upper_pm)) %>%
  data.table()

# Display minimum and maximum
min(data$pm25)
max(data$pm25)

# Make sure data is stored in correct form 
data[, region := as.factor(region)]
data[, race := as.factor(race)]
data[, entry_age_break := as.factor(entry_age_break)]
data[, dual := as.factor(dual)]
data[, year := factor(year)]
data[, followup_year := factor(followup_year)]

# histogram of mortality rate
data_hist <- 
  data %>%
  ggplot() + 
  geom_histogram(aes(x = mort_rate)) + 
  labs(x = "Mortaliry rate") + 
  theme_bw()

ggsave(data_hist, file = paste0(out_dir, "/hist_mortality.pdf"))


# Create histogram of PM2.5 concentration
gg_pm <- 
  data %>%
  ggplot() + 
  geom_histogram(aes(x = pm25)) + 
  labs(x = "PM concentration") + 
  theme_bw()

ggsave(gg_pm, file = paste0(out_dir, "/pm_concentration.pdf"))


# Replace zero with half of smallest value
min_positive_rate <- data[mort_rate > 0][mort_rate == min(mort_rate), mort_rate]
replacement_rate <- min_positive_rate / 2

# Now replace those with 0 
data[mort_rate == 0, mort_rate := replacement_rate]
data[, log_mort := log(mort_rate)]

# Create histogram of mortality rate
gg_log_mort <- 
  data %>%
  ggplot() + 
  geom_histogram(aes(x = log_mort)) + 
  labs(x = "Log Mortaliry rate") + 
  theme_bw()

ggsave(gg_log_mort, file = paste0(out_dir, "/log_mort.pdf"))


# Save input data to model folder (make sure not pushed to github)
saveRDS(data, file = paste0(out_dir, "/input_data.RDS"))

# Now create m out of n dataset for uncertainty quantification 
#(block is at the zip code level)
zip_codes <- unique(data$zip)
n <- length(zip_codes)
m <- round(n / log(n))

dir.create(paste0(out_dir, "/boostrap/"))

# Save 100 datasets for bootstrap
lapply(1:100, function(i) {
  set.seed(i)
  resample_zip_codes = sample(zip_codes, m, replace = T) # Sample zip with replacement
  resampled_data <- data[zip %in% resample_zip_codes] # save those zip codes
  saveRDS(resampled_data, file = paste0(out_dir, "/boostrap/boot_data_", i, ".RDS"))
})

message("Done")


