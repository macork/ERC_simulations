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

# # Load causal fit 
# causal_fit <- readRDS(file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/causal_fit.RDS"))
# causal_fit <- readRDS(file = paste0(proj_dir, "ERC_simulation/data_application/model_fits/causal_fit.RDS"))

# find min and max
min_pm <- min(data$pm25_ensemble)
max_pm <- max(data$pm25_ensemble)

data_prediction <- 
  rbindlist(lapply(seq(min_pm, max_pm, length.out = 100), function(pot_exp) {
    
    potential_data <- mutate(data, pm25_ensemble = pot_exp)
    
    # Fit each model and take mean for ERC
    potential_outcome <- 
      potential_data %>% 
      mutate(
        linear_model = predict(linear_fit, potential_data, type = "response"),
        gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
        #change_model = predict(change_model, newdata = potential_data, type = "response"),
        linear_ent = predict(linear_ent, potential_data, type = "response"),
        gam_ent = predict(gam_ent, newdata = potential_data, type = "response")
        #change_ent = predict(change_ent, newdata = potential_data, type = "response"),
        #causal_model =  predict(causal_gps, newdata = rename(potential_data, w = pm25_ensemble), type = "response")
             ) %>% 
      dplyr::select(pm25_ensemble, linear_model, linear_ent, gam_model, gam_ent) %>% 
      summarize_all(mean) %>% 
      data.table()
    return(potential_outcome)
  }))

saveRDS(data_prediction, file = paste0(model_dir, "data_prediction.RDS"))

model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent")

# Plot the fits from this sample
plot_fits <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(pm25_ensemble) %>%
  ggplot(aes(x = pm25_ensemble, y = exp(prediction), color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Relative risk of death")


ggsave(plot_fits, file = paste0(model_dir, "plot_fit.pdf"))
