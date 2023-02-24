# Script to run after models are fit, aggregates and creates curve
library(data.table)
library(tidyverse)
library(WeightIt)
library(fastDummies)
library(mgcv)
library(Rcpp)
library(RcppEigen)
library(CausalGPS)


# Get command arguments (first argument is which input data to use, second is name for model run)
args <- commandArgs(T)
data_flag <- as.character(args[[1]]) #supply tag for input data
model_flag = as.character(args[[2]])

# For publication tags are listed below
data_flag = "adjust_quantile"
model_flag = "old_boot"

# set directories
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/ERC_simulation"
model_dir <- paste0(proj_dir, "/Simulation_studies/data_application/model_fits/", model_flag, "/")

# Load in input data 
data <- readRDS(paste0(proj_dir, "/Medicare_data/model_input/", data_flag, "/input_data.RDS"))

#  If needed, load models 
# linear_fit <- readRDS(file = paste0(model_dir, "linear.RDS"))
# gam_fit <- readRDS(file = paste0(model_dir, "gam.RDS"))
# linear_ent <- readRDS(file = paste0(model_dir, "linear_ent.RDS"))
# gam_ent <- readRDS(file = paste0(model_dir, "gam_ent.RDS"))
# 
# # Load causal fit 
# pseudo_pop <- readRDS(file = paste0(model_dir, "pseudo_pop.RDS"))
# causal_gps <- readRDS(file = paste0(model_dir, "causal_fit.RDS"))

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

# Load data predictions
data_prediction <- readRDS(file = paste0(model_dir, "data_prediction.RDS"))

# Produce estimates of mortality rates
mortality_data <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25_ensemble) %>% 
  mutate(prediction = exp(prediction))

saveRDS(mortality_data, file = paste0(model_dir, "mortality_data.RDS"))

# Create plot of mean mortality rate for each model
mortality_plot <- 
  mortality_data %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  theme_bw() + 
  labs(x = "Exposure concentration", y = "Mortality rate")

ggsave(mortality_plot, file = paste0(model_dir, "mortality_plot.pdf"))

# Create relative rate by comparing mortality to underlying mortality at 12 ug/m3
twelve_index <- which.min(abs(data_prediction$pm25_ensemble - 12))

# Create dataset of relative rate
relative_rate_data <- 
  data_prediction %>% 
  arrange(pm25_ensemble) %>%
  mutate(linear_model = linear_model - data_prediction[twelve_index, linear_model]) %>%
  mutate(linear_ent = linear_ent - data_prediction[twelve_index, linear_ent]) %>%
  mutate(gam_model = gam_model - data_prediction[twelve_index, gam_model]) %>%
  mutate(gam_ent = gam_ent - data_prediction[twelve_index, gam_ent]) %>%
  mutate(causal_model = causal_model - data_prediction[twelve_index, causal_model]) %>%
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  mutate(prediction = exp(prediction)) %>% 
  data.table()

# Create plot of mean relative rate
relative_rate_plot <-
  relative_rate_data %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  coord_cartesian(xlim = c(min(data_prediction$pm25_ensemble), 12)) + 
  theme_bw() + 
  labs(x = "Annual Average PM2.5", y = "Relative rate of mortality")

ggsave(relative_rate_plot, file = paste0(model_dir, "relative_rate_plot.pdf"))

# Load contrast data (plots contrast of entire population at 12 ug/m3 to lower standard)
contrast_data <- readRDS(file = paste0(model_dir, "contrast_data.RDS"))


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
  group_by(pm25_ensemble, model) %>%
  summarize(var = var(prediction) * m/n) %>% # Scale variance, get se
  mutate(se = sqrt(var))

# Plot of mortality rate with standard errors
plot_label <- c("causal_model" = "CausalGPS", "gam_ent" = "GAM entropy", "gam_model" = "GAM",
                "linear_ent" = "Linear entropy", "linear_model" = "Linear")
gg_mort_se <- 
  mortality_data %>%
  left_join(mort_se, by = c("pm25_ensemble", "model")) %>%
  mutate(upper = prediction + 1.96*se, 
         lower = prediction - 1.96*se) %>%
  mutate(prediction = 100000 * prediction,
         lower = 100000 * lower, 
         upper = 100000 * upper) %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  geom_ribbon((aes(ymin = lower, ymax = upper, fill = model)), alpha = 0.4, color = NA) + 
  labs(x = "Annual Average PM2.5", y = "Mortality rate per\n100,000 person-years") + 
  theme_bw(base_size = 14) + 
  scale_fill_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_linetype_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_color_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  #coord_cartesian(xlim = c(4, 12), ylim = c(0.7, 1.2)) + 
  theme(legend.text=element_text(size=12))
  
ggsave(gg_mort_se, file = paste0(model_dir, "mortality_se.pdf"), width = 12, height = 8)

# Now calculate uncertainty for relative risk 
# scramble index so when dividing by risk at 12 we geg uncertainty interval at 12 as well 
index <- sample(1:100, replace = F)
pm_combaritor <- data_prediction[twelve_index, pm25_ensemble] # pm comparitor, should be very close to 12

# Get boot estimates of relative rate
relative_rate_boot <- 
  rbindlist(lapply(1:100, function(i){
    selected_sim <- mort_boot[sim == i] # grab simulation
    # Get reference mortality rate from another simulation
    reference_pm <- 
      mort_boot[sim == index[i] & pm25_ensemble == pm_combaritor, ] %>% 
      rename(reference = prediction) %>% select(model, reference)
    
    # Divide mortality rate by reference 
    relative_rate_output <- 
      selected_sim %>% 
      left_join(reference_pm, by = "model") %>% # Add reference
      mutate(relative_rate = prediction / reference) %>% 
      select(pm25_ensemble, model, relative_rate, sim)
    
    return(relative_rate_output)
  }))

# Now calculate standard error at each concentration
relative_rate_se <- 
  relative_rate_boot %>%
  group_by(pm25_ensemble, model) %>%
  summarize(var = var(relative_rate) * m/n) %>%
  mutate(se = sqrt(var))

gg_relative_se <- 
  relative_rate_data %>%
  left_join(relative_rate_se, by = c("pm25_ensemble", "model")) %>%
  mutate(upper = prediction + 1.96*se, 
         lower = prediction - 1.96*se) %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  geom_ribbon((aes(ymin = lower, ymax = upper, fill = model)), alpha = 0.4, color = NA) + 
  coord_cartesian(xlim = c(3, 12), ylim = c(0.7, 1.1)) + 
  labs(x = "Annual Average PM2.5", y = "Relative mortality rate") + 
  theme_bw(base_size = 14) + 
  scale_fill_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_linetype_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_color_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_x_continuous(breaks = c(4, 6, 8, 10, 12)) + 
  #coord_cartesian(xlim = c(4, 12), ylim = c(0.7, 1.2)) + 
  theme(legend.text=element_text(size=12))

ggsave(gg_relative_se, file = paste0(model_dir, "relative_se.pdf"), width = 12, height = 8)


# Now generate standard error of contrasts
contrast_boot <- 
  rbindlist(lapply(1:100, function(i){
    boot_data <- readRDS(paste0(boot_dir, "/", i, "/contrast_data.RDS"))
    boot_data <- data.table(boot_data)
    return(boot_data)
  }))

contrast_se <- 
  contrast_boot %>%
  pivot_longer(-contrast, "model", "value") %>%
  group_by(contrast, model) %>% 
  summarize(var = var(value) * m/n) %>%
  mutate(se = sqrt(var))

gg_contrast_se <- 
  contrast_data %>%
  pivot_longer(-contrast, "model", "value") %>%
  left_join(contrast_se) %>%
  mutate(upper = value + 1.96*se, 
         lower = value - 1.96*se) %>%
  ggplot(aes(x = contrast, y = 100 * value, color = model)) +
  #geom_point(size = 0.2, position = position_dodge(width = 0.7)) +
  geom_pointrange(aes(ymin = 100 * lower, ymax = 100 * upper, fill = model), 
                  size = 0.25, position = position_dodge(width = 0.15)) + 
  theme_bw() + 
  scale_fill_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_linetype_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  scale_color_discrete("", labels = plot_label, guide = guide_legend(reverse=TRUE)) + 
  labs(x = "Contrast", y = "Percent reduction in mortality")

ggsave(gg_contrast_se, file = paste0(model_dir, "contrast_se.pdf"), width = 12, height = 8)


# Create plots of covariate balance ------------------------------------------------------
# Load entropy weights from model directory
entropy_weights <- readRDS(file = paste0(model_dir, "entropy_weights.RDS"))

# Get covariate balance pre and post (can take)
post_weight <- cov.wt(entropy_weights[, 1:{ncol(entropy_weights)-1}], entropy_weights$ent_weights, TRUE)$cor[, 1]
pre_weight <- cor(entropy_weights[, 1:{ncol(entropy_weights)-1}])[, 1]

# Write function to get pre and post weights


post_weight <- post_weight[-1] # Remove pm
post_dt <- 
  data.table(value = post_weight) %>%
  mutate(variable = names(post_weight), cor = "post")

pre_weight <- pre_weight[-1] # Remove pm
pre_dt <- 
  data.table(value = pre_weight) %>%
  mutate(variable = names(pre_weight), cor = "pre") %>%
  arrange(value)

correlation_table <- rbind(pre_dt, post_dt)

# Combine factors for plot
region_breaks <- data.table("region", c("regionSOUTH", "regionNORTHEAST", "regionWEST", "regionSOUTH"))
year_breaks <- data.table("year", paste0("year", 2000:2016))
follow_breaks <- data.table("follow_up_year", paste0("followup_year", 1:18))
entry_breaks <- data.table("entry_age_break", paste0("entry_age_break", c(1, 3, 5, 7)))
race_breaks <- data.table("race", paste0("race", 1:5))
factor_index <- 
  rbind(data.table(c("education", "medianhousevalue", "medhouseholdincome", "poverty",
             "pct_owner_occ", "dual", "sex"), 
           c("education", "medianhousevalue", "medhouseholdincome", "poverty",
                                                "pct_owner_occ", "dual1", "sex1")),
      region_breaks, year_breaks, follow_breaks, entry_breaks, race_breaks)

names(factor_index) = c("var", "variable")

correlation_table_merged <- 
  correlation_table %>% 
  left_join(factor_index) %>% 
  group_by(var, cor) %>%
  summarize(value = mean(value))

gg_correlation_entropy <-
  correlation_table_merged %>% 
  mutate(value = abs(value)) %>%
  mutate(var = factor(var, levels = c("year", "education", "entry_age_break", 
                                      "region", "medianhousevalue", "medhouseholdincome", 
                                      "poverty", "pct_owner_occ", "dual", "race", "sex",
                                      "follow_up_year"))) %>%
  ggplot(aes(x = var, y = value, color = cor, group = cor)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_point() + 
  geom_line(size = 0.3) +
  theme_bw() +
  coord_flip() + 
  labs(y = "Absolute correlation", x = "Covariate", title = "Covariate balance entropy based weights")

ggsave(gg_correlation_entropy, file = paste0(model_dir, "correlation_entropy.pdf"))

# Now plot covariate balance for causalGPS
pseudo_pop <- readRDS(file = paste0(model_dir, "pseudo_pop.RDS"))
balance_gps <- plot(pseudo_pop)
ggsave(balance_gps, file = paste0(model_dir, "balance_gps.pdf"))
