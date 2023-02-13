# Code to prep data for model fitting (only need to run once for each new decision)
library(data.table)
library(tidyverse)
library(chngpt)
library(WeightIt)
library(fastDummies)
library(mgcv)
library(Rcpp)
library(RcppEigen)
library("CausalGPS", lib.loc = "/n/home_fasse/mcork/apps/ERC_simulation/R_4.0.5")

input_flag <- "adjust_quantile"

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
out_dir <- paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/")
dir.create(out_dir)

data <- readRDS(paste0(proj_dir, "analytic/aggregated_2000-2016_medicare_mortality_pm25_zip/aggregate_data.RDS"))

#data_test <- sample_n(data, 5000)
# Aggregate age groups and races 
# Only include, white, black, hispanic, other 
#data[, race := ifelse(race == 0, race = 1, race)]

# Aggreate age ranges to 10 year ranges
# 1 is 65-74
data[, entry_age_break := ifelse(entry_age_break == 2, 1, entry_age_break)] 
# 3 is 75-84
data[, entry_age_break := ifelse(entry_age_break == 4, 3, entry_age_break)]
# 5 is 85-94
data[, entry_age_break := ifelse(entry_age_break == 6, 5, entry_age_break)]
# 7 is 95 plus
data[, entry_age_break := ifelse(entry_age_break == 8, 7, entry_age_break)]

# See percent of age in each category
count(data, entry_age_break) %>% mutate(n = n / sum(n))

# Now aggregate over races (removing unknown race)
data <- 
  data %>% 
  filter(race != 0) %>% 
  mutate(race = ifelse(race == 6, 3, race))

# See percent of each race in category
count(data, race) %>% mutate(n = n / sum(n))

# now aggregate over age, sex, dual, zip code, year
group_cols <- setdiff(names(data), c("dead", "time_count"))

#data1 <- data %>% select(-dead, -time_count, -followup_year) %>% data.table()

# Collapse over these values in data table to save time
data <- data[, .(dead = sum(dead), time_count = sum(time_count)), by = group_cols]

# Find minimum 
# Now perform some functions of the data as a data exploration
#data <- aggregate_data[sample(.N, 1000)]
data[, mort_rate := dead / time_count]

# now trim upper and lower quantiles for data analysis
lower_pm <- quantile(data$pm25_ensemble, 0.025)
upper_pm <- quantile(data$pm25_ensemble, 0.975)

data <-
  data %>%
  filter(between(pm25_ensemble, lower_pm, upper_pm)) %>%
  data.table()

# # Create weighted percentile by the number of counts
# data1 <- 
#   data %>% 
#   arrange(pm25_ensemble) %>% 
#   mutate(wtd_ptile = lag(cumsum(time_count), default = 0)/(sum(time_count) - 1))
# 
# wtd_lower <- data1 %>% filter(wtd_ptile < .025) %>% slice_tail() %>% pull(pm25_ensemble)
# wtd_upper <- data1 %>% filter(wtd_ptile >= .975) %>% slice_head()

# Display minimum and maximum
min(data$pm25_ensemble)
max(data$pm25_ensemble)

# Make sure data is stored in correct form 
data[, region := as.factor(region)]
data[, sex := as.factor(sex - 1)]
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


# Now look at PM concentraion
#quantile(aggregate_data$pm25_ensemble)

gg_pm <- 
  data %>%
  ggplot() + 
  geom_histogram(aes(x = pm25_ensemble)) + 
  labs(x = "PM concentration") + 
  theme_bw()

ggsave(gg_pm, file = paste0(out_dir, "/pm_concentration.pdf"))


# Replace zero with half of smallest value
min_positive_rate <- data[mort_rate > 0][mort_rate == min(mort_rate), mort_rate]
replacement_rate <- min_positive_rate / 2

# Now replace those with 0 
data[mort_rate == 0, mort_rate := replacement_rate]
data[, log_mort := log(mort_rate)]

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
n <- nrow(data)
m <- round(n / log(n))

dir.create(paste0(out_dir, "/boostrap/"))

# Save 100 datasets for bootstrap, will later be used
lapply(1:100, function(i) {
  set.seed(i)
  resampled_data <- data[sample(.N, m)]
  saveRDS(resampled_data, file = paste0(out_dir, "/boostrap/boot_data_", i, ".RDS"))
})

message("Done")


