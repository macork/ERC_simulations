# Script to run after finished running through simulation

library(tidyverse)
library(data.table)

exp_relationship = "linear"
adjust_confounder = F
time_stamp = "0803_14"
replicates = 100

results_dir <- paste0("/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/results/",
                      exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")



metrics <- 
  rbindlist(lapply(1:replicates, function(rep){
    sim_metrics <- try(readRDS(paste0(results_dir, "sim_final_", rep, ".RDS"))[[1]], silent = T)
    if (is.data.frame(sim_metrics)) {
      data.table(sim_metrics)
      return(sim_metrics)
    } else {
      return()
    }
  }))

predictions <- 
  rbindlist(lapply(1:replicates, function(rep){
    sim_pred <- try(readRDS(paste0(results_dir, "sim_final_", rep, ".RDS"))[[2]], silent = T)
    if (is.data.frame(sim_pred)) {
      data.table(sim_pred)
      return(sim_pred)
    } else {
      return()
    }
  }))

correlation <- 
  rbindlist(lapply(1:replicates, function(rep){
    sim_corr <- try(readRDS(paste0(results_dir, "sim_final_", rep, ".RDS"))[[3]], silent = T)
    if (is.data.frame(sim_corr)) {
      data.table(sim_corr)
      return(sim_corr)
    } else {
      return()
    }
  }))

saveRDS(list(metrics = metrics, predictions = predictions, correlation = correlation),
        file = paste0(results_dir, "aggregated_results.RDS"))


gg_correlation_ipw <-
  correlation %>% 
  filter(method == "ipw") %>% 
  group_by(gps_mod, covariate) %>% 
  dplyr::summarize(pre_weight = mean(pre_cor), post_weight = mean(post_cor)) %>% 
  pivot_longer(c("pre_weight", "post_weight")) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  ggplot(aes(x = covariate, y = value, color = name, group = name)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  coord_flip() + 
  facet_wrap(~ gps_mod) + 
  labs(y = "Absolute correlation", x = "Covariate", title = "Comparing covariate balance under different design settings for IPW")
gg_correlation_ipw

gg_correlation_default <-
  correlation %>% 
  filter(method == "causal_gps_default") %>% 
  group_by(gps_mod, covariate) %>% 
  dplyr::summarize(pre_weight = mean(pre_cor), post_weight = mean(post_cor)) %>% 
  pivot_longer(c("pre_weight", "post_weight")) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  ggplot(aes(x = covariate, y = value, color = name, group = name)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  coord_flip() + 
  facet_wrap(~ gps_mod) + 
  labs(y = "Absolute correlation", x = "Covariate", title = "Comparing covariate balance under different design settings for causal GPS default")

gg_correlation_default

gg_correlation_tuned <-
  correlation %>% 
  filter(method == "causal_gps_tuned") %>% 
  group_by(gps_mod, covariate) %>% 
  dplyr::summarize(pre_weight = mean(pre_cor), post_weight = mean(post_cor)) %>% 
  pivot_longer(c("pre_weight", "post_weight")) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  ggplot(aes(x = covariate, y = value, color = name, group = name)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  coord_flip() + 
  facet_wrap(~ gps_mod) + 
  labs(y = "Absolute correlation", x = "Covariate", title = "Comparing covariate balance under different design settings for causal GPS default")

gg_correlation_tuned

plot_data <- 
  metrics %>% 
  mutate(absolute_bias = abs(bias)) %>% 
  tidyr::pivot_longer(c("bias", "absolute_bias", "mse")) %>% 
  mutate(model = factor(model, levels = c("linear_model", "linear_gps", "gam_model", "gam_gps", "causal_gps_default", "causal_gps_tuned", "change_model", "propensity_change")))  %>% 
  dplyr::group_by(model, gps_mod, sample_size, name) %>% 
  dplyr::summarize(mean = mean(value), 
                   lower = quantile(value, 0.05),
                   upper = quantile(value, 0.95)) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction")))

plot_data %>% 
  filter(name == "bias") %>% 
  ggplot() + 
  geom_pointrange(aes(y = mean, ymin = lower, ymax = upper, x = factor(gps_mod), color = model), position = position_dodge(.6)) +
  geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
  labs(x = "Scenario", y = "Value") + 
  theme_bw() + 
  coord_cartesian(ylim = c(-20, 20)) + 
  facet_wrap(sample_size ~ name, scales = "free")

plot_data %>% 
  filter(name == "absolute_bias", sample_size == 1000) %>% 
  ggplot() + 
  geom_pointrange(aes(y = mean, ymin = lower, ymax = upper, x = factor(gps_mod), color = model), position = position_dodge(.6)) +
  geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
  labs(x = "Scenario", y = "Value") + 
  coord_cartesian(ylim = c(0, 35)) + 
  theme_bw() + 
  facet_grid(sample_size ~ name, scales = "free")

plot_data %>% 
  filter(name == "mse") %>% 
  ggplot() + 
  geom_pointrange(aes(y = mean, ymin = lower, ymax = upper, x = factor(gps_mod), color = model), position = position_dodge(.6)) +
  geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype = "dashed", size = 0.2) + 
  labs(x = "Scenario", y = "Value") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 400)) + 
  facet_grid(sample_size ~ name, scales = "free")

