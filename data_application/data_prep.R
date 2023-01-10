# Code to prep data for model fitting (only need to run once for each new decision)
library(data.table)
library(tidyverse)
library(chngpt)
library(WeightIt)
library(fastDummies)
library(mgcv)
library(Rcpp)
library(RcppEigen)
library(CausalGPS)

input_flag <- "quantile_removed"

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
out_dir <- paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/")
dir.create(out_dir)

data <- readRDS(paste0(proj_dir, "analytic/aggregated_2000-2016_medicare_mortality_pm25_zip/aggregate_data.RDS"))

#data_test <- sample_n(data, 1000)
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

# Now aggregate over races
data <- 
  data %>% 
  filter(race != 0) %>% 
  mutate(race = ifelse(race == 6, 3, race))

# See percent of each race in category
count(data, race) %>% mutate(n = n / sum(n))

# now aggregate over age, sex, dual, zip code, year
group_cols <- setdiff(names(data), c("dead", "time_count"))

# Collapse over these values in data table to save time
data <- data[, .(dead = sum(dead), time_count = sum(time_count)), by = group_cols]

# Find minimum 
# Now perform some functions of the data as a data exploration
#data <- aggregate_data[sample(.N, 1000)]
data[, mort_rate := dead / time_count]

# now trim upper and lower quantiles
lower_pm <- quantile(data$pm25_ensemble, 0.025)
upper_pm <- quantile(data$pm25_ensemble, 0.975)

data <- 
  data %>% 
  filter(between(pm25_ensemble, lower_pm, upper_pm)) %>% 
  data.table()

# Find minimum and maximum
min_pm <- min(data$pm25_ensemble)
max_pm <- max(data$pm25_ensemble)

# Make sure data is stored in correct form 
data[, region := as.factor(region)]
data[, sex := as.factor(sex)]
data[, race := as.factor(race)]
data[, entry_age_break := as.factor(entry_age_break)]
data[, dual := as.factor(dual)]
names(data)

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
