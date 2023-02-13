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
input_flag <- as.character(args[[1]])
model_flag = as.character(args[[2]])

# load in data
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/"
data <- readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", input_flag, "/input_data.RDS"))
model_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")

# load the models if doing post processing only
linear_fit <- readRDS(file = paste0(model_dir, "linear.RDS"))
gam_fit <- readRDS(file = paste0(model_dir, "gam.RDS"))
linear_ent <- readRDS(file = paste0(model_dir, "linear_ent.RDS"))
gam_ent <- readRDS(file = paste0(model_dir, "gam_ent.RDS"))

# Load causal fit 
pseudo_pop <- readRDS(file = paste0(model_dir, "pseudo_pop.RDS"))
causal_gps <- readRDS(file = paste0(model_dir, "causal_fit.RDS"))

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

data_prediction <- readRDS(file = paste0(model_dir, "data_prediction.RDS"))

mortality_data <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25_ensemble) %>% 
  mutate(prediction = exp(prediction))

saveRDS(mortality_data, file = paste0(model_dir, "mortality_data.RDS"))

mortality_plot <- 
  mortality_data %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Mortality rate")

ggsave(mortality_plot, file = paste0(model_dir, "mortality_plot.pdf"))


relative_rate_data <- 
  data_prediction %>% 
  arrange(pm25_ensemble) %>%
  mutate(linear_model = linear_model - data_prediction$linear_model[1]) %>%
  mutate(linear_ent = linear_ent - data_prediction$linear_ent[1]) %>%
  mutate(gam_model = gam_model - data_prediction$gam_model[1]) %>%
  mutate(gam_ent = gam_ent - data_prediction$gam_ent[1]) %>%
  mutate(causal_model = causal_model - data_prediction$causal_model[1]) %>%
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  mutate(prediction = exp(prediction)) %>% 
  data.table()

saveRDS(relative_rate_data, file = paste0(model_dir, "relative_rate.RDS"))

relative_rate_data2 <- 
  data_prediction %>% 
  arrange(pm25_ensemble) %>%
  mutate(linear_model = linear_model - data_prediction[70, linear_model]) %>%
  mutate(linear_ent = linear_ent - data_prediction[70, linear_ent]) %>%
  mutate(gam_model = gam_model - data_prediction[70, gam_model]) %>%
  mutate(gam_ent = gam_ent - data_prediction[70,gam_ent]) %>%
  mutate(causal_model = causal_model - data_prediction[70,causal_model]) %>%
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  mutate(prediction = exp(prediction)) %>% 
  data.table()

relative_rate_plot <-
  relative_rate_data2 %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  coord_cartesian(xlim = c(4, 12)) + 
  labs(x = "Annual Average PM2.5", y = "Hazard Ratio")

ggsave(relative_rate_plot, file = paste0(model_dir, "relative_rate_plot.pdf"))

# Load contrast data
contrast_data <- readRDS(file = paste0(model_dir, "contrast_data.RDS"))


# Load boot data -------------------------------------------------------------

m <- nrow(readRDS(paste0(proj_dir, "ERC_simulation/Medicare_data/model_input/", 
                         input_flag, "/boostrap/boot_data_100.RDS")))
n <- nrow(data)

boot_dir <- paste0(model_dir, "/boostrap/")
mort_boot <- 
  rbindlist(lapply(1:100, function(i){
    boot_data <- readRDS(paste0(boot_dir, "/", i, "/mortality_data.RDS"))
    boot_data <- data.table(boot_data)
    boot_data[, sim := i]
    return(boot_data)
  }))

mort_se <- 
  mort_boot %>%
  group_by(pm25_ensemble, model) %>%
  summarize(var = var(prediction) * m/n) %>%
  mutate(se = sqrt(var))

gg_mort_se <- 
  mortality_data %>%
  left_join(mort_se) %>%
  mutate(upper = prediction + 1.96*se, 
         lower = prediction - 1.96*se) %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  geom_ribbon((aes(ymin = lower, ymax = upper, fill = model)), alpha = 0.4, color = NA) + 
  labs(x = "Exposure concentration", y = "Mortality rate")
  
ggsave(gg_mort_se, file = paste0(model_dir, "mortality_se.pdf"))


#mortality_data$upper <- mortality_data$prediction + 1.96 * mort_se$se

relative_boot <- 
  rbindlist(lapply(1:100, function(i){
    boot_data <- readRDS(paste0(boot_dir, "/", i, "/relative_rate.RDS"))
    boot_data <- data.table(boot_data)
    boot_data[, sim := i]
    return(boot_data)
  }))

relative_rate_se <- 
  relative_boot %>%
  group_by(pm25_ensemble, model) %>%
  summarize(var = var(prediction) * m/n) %>%
  mutate(se = sqrt(var))


gg_relative_se <- 
  relative_rate_data2 %>%
  left_join(relative_rate_se) %>%
  mutate(upper = prediction + 1.96*se, 
         lower = prediction - 1.96*se) %>%
  ggplot(aes(x = pm25_ensemble, y = prediction, color = model, linetype = model)) +
  geom_line() +
  geom_ribbon((aes(ymin = lower, ymax = upper, fill = model)), alpha = 0.4, color = NA) + 
  theme_bw() + 
  coord_cartesian(xlim = c(4, 12)) + 
  labs(x = "Annual Average PM2.5", y = "Hazard Ratio")

ggsave(gg_relative_se, file = paste0(model_dir, "relative_se.pdf"))

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
  labs(x = "Contrast", y = "Percent reductions")


# Now work on covariate balance
entropy_weights <- readRDS(file = paste0(model_dir, "entropy_weights.RDS"))
post_weight <- cov.wt(entropy_weights[, 1:20], entropy_weights$ent_weights, TRUE)$cor[, 1]
pre_weight <- cor(entropy_weights[, 1:20])[, 1]

post_weight <- post_weight[-1]
post_dt <- 
  data.table(value = post_weight) %>%
  mutate(variable = names(post_weight), cor = "post")

pre_weight <- pre_weight[-1]
pre_dt <- 
  data.table(value = pre_weight) %>%
  mutate(variable = names(pre_weight), cor = "pre") %>%
  arrange(value)

correlation_table <- rbind(pre_dt, post_dt)

data.table(post_weight[-1]) %>% pivot_longer()


gg_correlation_entropy <-
  correlation_table %>% 
  mutate(value = abs(value)) %>%
  #arrange(value) %>%
  ggplot(aes(x = variable, y = value, color = cor, group = cor)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_point() + 
  geom_line(size = 0.3) +
  theme_bw() +
  coord_flip() + 
  labs(y = "Absolute correlation", x = "Covariate", title = "Covariate balance entropy based weights")

ggsave(gg_correlation_entropy, file = paste0(model_dir, "correlation_entropy.pdf"))

# Now plot covariate balance
balance_gps <- plot(pseudo_pop)
ggsave(balance_gps, file = paste0(model_dir, "balance_gps.pdf"))
