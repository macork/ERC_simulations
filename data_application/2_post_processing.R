# Script to run after models are fit, aggregates and creates plots
rm(list = ls())

# Load libraries
library(data.table)
library(tidyverse)
library(WeightIt)
library(fastDummies)
library(mgcv)
library(Rcpp)
library(RcppEigen)
library(CausalGPS)
library(cobalt)


# Get command arguments (first argument is which input data to use, second is name for model run)
args <- commandArgs(T)
data_flag <- as.character(args[[1]]) #supply tag for input data
model_flag = as.character(args[[2]])

# For publication tags are listed below
data_flag = "kevin_trim_90"
model_flag = "Causal_by_year"

# set directories
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/ERC_simulation"
model_dir <- paste0(proj_dir, "/Simulation_studies/data_application/model_fits/", model_flag, "/")

# Load in input data 
data <- readRDS(paste0(proj_dir, "/Medicare_data/model_input/", data_flag, "/input_data.RDS"))

#  If needed, load models
linear_fit <- readRDS(file = paste0(model_dir, "linear.RDS"))
gam_fit <- readRDS(file = paste0(model_dir, "gam.RDS"))
linear_ent <- readRDS(file = paste0(model_dir, "linear_ent.RDS"))
gam_ent <- readRDS(file = paste0(model_dir, "gam_ent.RDS"))

# Load causal fit
pseudo_pop <- readRDS(file = paste0(model_dir, "pseudo_pop.RDS"))
causal_gps <- readRDS(file = paste0(model_dir, "causal_fit.RDS"))

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

# Load data predictions
data_prediction <- readRDS(file = paste0(model_dir, "data_prediction.RDS"))

# Produce estimates of mortality rates
mortality_data <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25)

saveRDS(mortality_data, file = paste0(model_dir, "mortality_data.RDS"))

# Create plot of mean mortality rate for each model
mortality_plot <- 
  mortality_data %>%
  ggplot(aes(x = pm25, y = exp(prediction), color = model, linetype = model)) +
  geom_line() +
  theme_bw() + 
  labs(x = "Exposure concentration", y = "Mortality rate")

ggsave(mortality_plot, file = paste0(model_dir, "mortality_plot.pdf"))

# Create relative rate by comparing mortality to underlying mortality at 12 ug/m3
twelve_index <- which.min(abs(data_prediction$pm25 - 12))

# Create dataset of relative rate
relative_rate_data <- 
  data_prediction %>% 
  arrange(pm25) %>%
  mutate(linear_model = linear_model - data_prediction[twelve_index, linear_model]) %>%
  mutate(linear_ent = linear_ent - data_prediction[twelve_index, linear_ent]) %>%
  mutate(gam_model = gam_model - data_prediction[twelve_index, gam_model]) %>%
  mutate(gam_ent = gam_ent - data_prediction[twelve_index, gam_ent]) %>%
  mutate(causal_model = causal_model - data_prediction[twelve_index, causal_model]) %>%
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  data.table()

# Create plot of mean relative rate
relative_rate_plot <-
  relative_rate_data %>%
  ggplot(aes(x = pm25, y = exp(prediction), color = model, linetype = model)) +
  geom_line() +
  coord_cartesian(xlim = c(min(data_prediction$pm25), 12)) + 
  theme_bw() + 
  labs(x = "Annual Average PM2.5", y = "Relative rate of mortality")

ggsave(relative_rate_plot, file = paste0(model_dir, "relative_rate_plot.pdf"))

# Load contrast data (plots contrast of entire population at 12 ug/m3 to lower standard)
#contrast_data <- readRDS(file = paste0(model_dir, "contrast_data.RDS"))

# Load boot data for uncertainty quantification --------------------------------------------------
# See size of m out of n bootstrap
m <- nrow(readRDS(paste0(proj_dir, "/Medicare_data/model_input/", data_flag, "/boostrap/boot_data_99.RDS")))
n <- nrow(data)

# Lood boostrap of mortality rates
boot_dir <- paste0(model_dir, "/boostrap/")
mort_boot <- 
  rbindlist(lapply(1:100, function(i){
    boot_data <- readRDS(paste0(boot_dir, "/", i, "/mortality_data.RDS"))
    boot_data <- data.table(boot_data)
    boot_data[, sim := i]
    return(boot_data)
  }))

# Calculate standard error at each concentration, multiply by m/n 
mort_se <- 
  mort_boot %>%
  group_by(pm25, model) %>%
  mutate(prediction = log(prediction)) %>% # Put back into log space
  summarize(var = var(prediction) * m/n) %>% # Scale variance, get se
  mutate(se = sqrt(var))

# Plot of mortality rate with standard errors
plot_label <- c("causal_model" = "CausalGPS", "gam_ent" = "GAM entropy", "gam_model" = "GAM",
                "linear_ent" = "Linear entropy", "linear_model" = "Linear")
gg_mort_se <- 
  mortality_data %>%
  left_join(mort_se, by = c("pm25", "model")) %>%
  mutate(upper = exp(prediction + 1.96*se), 
         lower = exp(prediction - 1.96*se)) %>%
  mutate(prediction = 100000 * exp(prediction),
         lower = 100000 * lower, 
         upper = 100000 * upper) %>% 
  ggplot(aes(x = pm25, y = prediction, color = model, linetype = model)) +
  geom_line() +
  geom_ribbon((aes(ymin = lower, ymax = upper, fill = model, linetype = model, color = model)), alpha = 0.2) + 
  #geom_ribbon((aes(ymin = lower, ymax = upper, fill = model)), alpha = 0.3, color = NA) + 
  labs(x = "Annual Average PM2.5", y = "Mortality rate per 100,000 person-years") + 
  theme_bw(base_size = 16) + 
  scale_fill_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_linetype_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_color_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  theme(legend.text=element_text(size=16),
        legend.position=c(.85,.25))
  #coord_cartesian(xlim = c(4, 12), ylim = c(0.7, 1.2)) + 
  
ggsave(gg_mort_se, 
       file = paste0(model_dir, "mortality_se2.png"), 
       width = 10, height = 7,
       dpi = 400)

# Now calculate uncertainty for relative risk 
# scramble index so when dividing by risk at 12 we get uncertainty interval at 12 as well 
index <- sample(1:100, replace = F)
pm_combaritor <- data_prediction[twelve_index, pm25] # pm comparitor, should be very close to 12

# Get boot estimates of relative rate
relative_rate_boot <- 
  rbindlist(lapply(1:100, function(i){
    selected_sim <- mort_boot[sim == i, ] %>% mutate(prediction = log(prediction)) 
    # Get reference mortality rate from another simulation
    reference_pm <- 
      mort_boot[sim == index[i] & pm25 == pm_combaritor, ] %>% 
      mutate(prediction = log(prediction)) %>% 
      rename(reference = prediction) %>% select(model, reference)
    
    # Divide mortality rate by reference 
    relative_rate_output <- 
      selected_sim %>% 
      left_join(reference_pm, by = "model") %>% # Add reference
      mutate(relative_rate = prediction - reference) %>% 
      select(pm25, model, relative_rate, sim)
    
    return(relative_rate_output)
  }))

# Now calculate standard error at each concentration
relative_rate_se <- 
  relative_rate_boot %>%
  group_by(pm25, model) %>%
  summarize(var = var(relative_rate) * m/n) %>%
  mutate(se = sqrt(var))

# plot of gg_relative_se
gg_relative_se <- 
  relative_rate_data %>%
  left_join(relative_rate_se, by = c("pm25", "model")) %>%
  mutate(upper = exp(prediction + 1.96*se), 
         lower = exp(prediction - 1.96*se)) %>%
  mutate(prediction = exp(prediction)) %>%
  ggplot(aes(x = pm25, y = prediction, color = model, linetype = model)) +
  geom_line() +
  geom_ribbon((aes(ymin = lower, ymax = upper, fill = model, linetype = model, color = model)), alpha = 0.2) + 
  coord_cartesian(xlim = c(6, 12), ylim = c(0.8, 1.05)) + 
  labs(x = "Annual Average PM2.5", y = "Relative mortality rate") + 
  theme_bw(base_size = 16) + 
  scale_fill_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_linetype_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_color_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) +
  theme(legend.text=element_text(size=16),
        legend.position=c(.85,.25))

ggsave(gg_relative_se, 
       file = paste0(model_dir, "relative_se_present.png"), 
       width = 10, height = 7,
       dpi = 400)


# Create plots of covariate balance ------------------------------------------------------
# Load entropy weights from model directory
# Assign variables as strata or confounders
strata_var <- c("female", "race", "dual", "entry_age_break", "followup_year")
confounders <- c("smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue", "mean_bmi",
                 "poverty", "education", "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax",
                 "winter_rmax", "year", "region")
entropy_data <- readRDS(file = paste0(model_dir, "entropy_weights.RDS"))
pseudo_pop <- readRDS(file = paste0(model_dir, "pseudo_pop.RDS"))


# Arrange entropy and pseudo population, should be same dataset with only different weights
# This is to make sure they are the same
arranged_psuedo <- pseudo_pop %>% arrange(Y, year) %>% rename(pm25 = w) %>% 
  select(!!! c(confounders, "pm25", "counter_weight")) %>% data.table()
arranged_entropy <- entropy_data %>% arrange(zip, year) %>%  select(!!! c(confounders, "pm25", "ent_weight")) %>% data.table()

# Test these data sets are equal other than weighting scheme
all.equal(select(arranged_psuedo, -counter_weight), select(arranged_entropy, -ent_weight))

# Grab correlation unbalance (nice wrapper for this)
unweighted_balance <- 
  bal.tab(reformulate(confounders, response = "pm25"),
          data = arranged_entropy,
          weights = arranged_entropy[["ent_weight"]],
          method = "weighting",
          stats = c("cor"),
          un = T,
          continuous = "std",
          #s.d.denom = "weighted",
          abs = T,
          thresholds = c(cor = .1), poly = 1)

unweighted_df <- 
  cbind(data.table(unweighted_balance$Balance), variable =  rownames(unweighted_balance$Balance)) %>% 
  select(variable, corr = Corr.Un) %>% mutate(type = "Unadjusted")

# Get weighted correlation for entropy balancing
ent_balance <- 
  bal.tab(reformulate(confounders, response = "pm25"),
          data = arranged_entropy,
          weights = arranged_entropy[["ent_weight"]],
          method = "weighting",
          stats = c("cor"),
          un = T,
          continuous = "std",
          s.d.denom = "weighted",
          abs = T,
          thresholds = c(cor = .1), poly = 1)

ent_df <- 
  cbind(data.table(ent_balance$Balance), variable =  rownames(ent_balance$Balance)) %>% 
  select(variable, corr = Corr.Adj) %>% mutate(type = "Entropy")

# Get weighted correlation for causalGPS
causal_balance <- 
  bal.tab(reformulate(confounders, response = "pm25"),
          data = arranged_psuedo,
          weights = arranged_psuedo[["counter_weight"]],
          method = "weighting",
          stats = c("cor"),
          un = T,
          continuous = "std",
          s.d.denom = "weighted",
          abs = T,
          thresholds = c(cor = .1), poly = 1)

causal_df <- 
  cbind(data.table(causal_balance$Balance), variable =  rownames(causal_balance$Balance)) %>% 
  select(variable, corr = Corr.Adj) %>% mutate(type = "CausalGPS")

# Bind together for balance table
balance_table <- rbind(unweighted_df, ent_df, causal_df) 

# Plot full balance table (not combining year and region)
balance_table %>%
  ggplot(aes(x = variable, y = corr, color = type, group = type)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_point() + 
  geom_line(linewidth = 0.3) +
  theme_bw(base_size = 14) + 
  theme(legend.text=element_text(size=12)) +
  coord_flip() + 
  scale_color_discrete("Correlation") + 
  labs(y = "Absolute correlation", x = "Covariates", title = "Covariate Balance Test")

# Now combine year and regions
region_breaks <- data.table(var = "region", 
                            variable = c("region_SOUTH", "region_NORTHEAST",
                                         "region_WEST", "region_MIDWEST"))
year_breaks <- data.table(var = "year", variable = paste0("year_", 2000:2016))

# Create way to merge all factors together
merge_factors <- 
  rbind(data.table(var = confounders, variable = confounders),
        region_breaks, year_breaks)

balance_merged <- 
  balance_table %>% 
  left_join(merge_factors) %>% 
  group_by(var, type) %>%
  summarize(corr = mean(corr))

# Order by decreasing unadjusted
order_var <- 
  balance_merged %>%
  filter(type == "Unadjusted") %>% 
  arrange(corr) %>% 
  pull(var)

mean_abs_corr <- 
  balance_table %>%
  group_by(type) %>%
  summarize(abs_cor = mean(corr)) %>% 
  data.table()
  
# Label covariates
cov_label = rev(c("% Below High School Education",
              "Summer Humidity",
              "% Black",
              "Region", 
              "Median Home Value", 
              "Summer Temperature",
              "Year",
              "Median Household Income",
              "Population Density",
              "% Below Poverty Level",
              "% Owner-occupied Housing",
              "Winter Temperature", 
              "% Ever Smoked",
              "% Hispanic", 
              "Winter Humidty", 
              "Mean BMI"))

# Create balance plot
balance_merged_plot <- 
  balance_merged %>% 
  mutate(var = factor(var, levels = order_var, labels = cov_label)) %>% 
  mutate(type = factor(type, levels = c("Unadjusted", "CausalGPS", "Entropy"))) %>%
  ggplot(aes(x = var, y = corr, color = type, group = type)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_hline(aes(yintercept = mean_abs_corr[type == "Unadjusted", abs_cor]), 
             linetype = "dotted", color = "red") + 
  geom_hline(aes(yintercept = mean_abs_corr[type == "CausalGPS", abs_cor]), 
             linetype = "dotted", color = "green") + 
  geom_hline(aes(yintercept = mean_abs_corr[type == "Entropy", abs_cor]), 
             linetype = "dotted", color = "blue") + 
  geom_point() + 
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  geom_line(linewidth = 0.3) +
  theme_bw(base_size = 16) + 
  theme(legend.text=element_text(size=12)) +
  theme(
    legend.position = c(0.835, 0.2),
    legend.title = element_blank(),
    legend.text=element_text(size=16)
  ) + 
  coord_flip() + 
  labs(y = "Absolute correlation", x = "Covariates", title = "")

ggsave(balance_merged_plot, 
       file = paste0(model_dir, "balance_plot.png"), width = 10, height = 10,
       dpi = 400)

# # Plots of contrasts (omit for now)
# contrast_boot <- 
#   rbindlist(lapply(1:100, function(i){
#     boot_data <- readRDS(paste0(boot_dir, "/", i, "/contrast_data.RDS"))
#     boot_data <- data.table(boot_data)
#     return(boot_data)
#   }))
# 
# contrast_se <- 
#   contrast_boot %>%
#   pivot_longer(-contrast, "model", "value") %>%
#   group_by(contrast, model) %>% 
#   summarize(var = var(value) * m/n) %>%
#   mutate(se = sqrt(var))
# 
# gg_contrast_se <- 
#   contrast_data %>%
#   pivot_longer(-contrast, "model", "value") %>%
#   left_join(contrast_se) %>%
#   mutate(upper = value + 1.96*se, 
#          lower = value - 1.96*se) %>%
#   ggplot(aes(x = contrast, y = 100 * value, color = model)) +
#   #geom_point(size = 0.2, position = position_dodge(width = 0.7)) +
#   geom_pointrange(aes(ymin = 100 * lower, ymax = 100 * upper, fill = model), 
#                   size = 0.25, position = position_dodge(width = 0.15)) + 
#   theme_bw() + 
#   scale_fill_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
#   scale_linetype_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
#   scale_color_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
#   labs(x = "Contrast", y = "Percent reduction in mortality")
# 
# ggsave(gg_contrast_se, file = paste0(model_dir, "contrast_se.pdf"), width = 12, height = 8)
