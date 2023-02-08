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

# load the linear model fit
linear_fit <- readRDS(file = paste0(model_dir, "linear.RDS"))

# Load GAM fit 
gam_fit <- readRDS(file = paste0(model_dir, "gam.RDS"))

# Load linear entropy weights
linear_ent <- readRDS(file = paste0(model_dir, "linear_ent.RDS"))

# Load gam entropy weights
gam_ent <- readRDS(file = paste0(model_dir, "gam_ent.RDS"))

# Load causal fit 
model_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/causal_run1/")
pseudo_pop <- readRDS(file = paste0(model_dir, "pseudo_pop5.RDS"))
psuedo_pop_frame <- data.table(pseudo_pop$pseudo_pop)

causal_gps <- 
  mgcv::bam(formula = Y ~ s(w, bs = 'cr', k = 4),
            family = "gaussian",
            data = psuedo_pop_frame,
            weights = counter_weight)


# causal_fit <- readRDS(file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/causal_fit.RDS"))
# causal_fit <- readRDS(file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/causal_fit.RDS"))

# find min and max
min_pm <- min(data$pm25_ensemble)
max_pm <- max(data$pm25_ensemble)

data_prediction <- 
  rbindlist(lapply(seq(min_pm, max_pm, length.out = 10), function(pot_exp) {
    
    potential_data <- mutate(data, pm25_ensemble = pot_exp)
    
    # Fit each model and take mean for ERC
    potential_outcome <- 
      potential_data %>% 
      mutate(
        linear_model = predict(linear_fit, potential_data, type = "response"),
        gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
        #change_model = predict(change_model, newdata = potential_data, type = "response"),
        linear_ent = predict(linear_ent, potential_data, type = "response"),
        gam_ent = predict(gam_ent, newdata = potential_data, type = "response"),
        #change_ent = predict(change_ent, newdata = potential_data, type = "response"),
        causal_model =  predict(causal_gps, newdata = rename(potential_data, w = pm25_ensemble), type = "response")
             ) %>% 
      dplyr::select(pm25_ensemble, linear_model, linear_ent, gam_model, gam_ent, causal_model) %>% 
      summarize_all(mean) %>% 
      data.table()
    return(potential_outcome)
  }))

saveRDS(data_prediction, file = paste0(model_dir, "data_prediction.RDS"))

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "causal_model")

data_prediction <- readRDS(file = paste0(model_dir, "data_prediction.RDS"))

# Plot the fits from this sample
mortality_plot <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25_ensemble) %>%
  ggplot(aes(x = pm25_ensemble, y = exp(prediction), color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Mortality rate")

ggsave(mortality_plot, file = paste0(model_dir, "mortality_plot.pdf"))

relative_rate_plot <- 
  data_prediction %>% 
  mutate(linear_model = linear_model - data_prediction$linear_model[1]) %>%
  mutate(linear_ent = linear_ent - data_prediction$linear_ent[1]) %>%
  mutate(gam_model = gam_model - data_prediction$gam_model[1]) %>%
  mutate(gam_ent = gam_ent - data_prediction$gam_ent[1]) %>%
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25_ensemble) %>%
  ggplot(aes(x = pm25_ensemble, y = exp(prediction), color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Relative mortality rate")

ggsave(relative_rate_plot, file = paste0(model_dir, "relative_rate_plot.pdf"))

# Plot the contrasts
contrasts <- c(8, 9, 10, 12)

b = c(8, 12)
c1 <- 
  data_prediction %>%  
  filter(abs(pm25_ensemble - b[1]) == min(abs(pm25_ensemble - b[1])))
c2 <- 
  data_prediction %>%  
  filter(abs(pm25_ensemble - b[2]) == min(abs(pm25_ensemble - b[2])))

contrast8_12 <- 1 - exp(c1 - c2) %>% select(!!model_types)

b = c(9, 12)
c1 <- 
  data_prediction %>%  
  filter(abs(pm25_ensemble - b[1]) == min(abs(pm25_ensemble - b[1])))
c2 <- 
  data_prediction %>%  
  filter(abs(pm25_ensemble - b[2]) == min(abs(pm25_ensemble - b[2])))

contrast9_12 <- 1 - exp(c1 - c2) %>% select(!!model_types)


b = c(10, 12)
c1 <- 
  data_prediction %>%  
  filter(abs(pm25_ensemble - b[1]) == min(abs(pm25_ensemble - b[1])))
c2 <- 
  data_prediction %>%  
  filter(abs(pm25_ensemble - b[2]) == min(abs(pm25_ensemble - b[2])))

contrast10_12 <- 1 - exp(c1 - c2) %>% select(!!model_types)

gg_contrast <- 
  rbind(contrast8_12, contrast9_12, contrast10_12) %>% 
  cbind(contrast = c(8, 9, 10)) %>% 
  data.table() %>% 
  pivot_longer(model_types) %>%
  ggplot() + 
  geom_point(aes(x = contrast, y = 100 * value, color = name), size = 3) +
  labs(x = "Lower emission limit", y = "Percent Reduction\nin relative mortality", 
       title = "Reduction in mortality by reducing limit from 12 ug/ml3")

ggsave(gg_contrast, file = paste0(model_dir, "contrast_plot.pdf"))

# Now work on covariate balance
model_dir <- paste0(proj_dir, "ERC_simulation/Simulation_studies/data_application/model_fits/", model_flag, "/")
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
