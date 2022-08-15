## File for creating figures from simulation results to present

library(tidyverse)
library(data.table)


exp_relationship = "sublinear"
adjust_confounder = T
time_stamp = "sublinear_interaction_complex_large"
replicates = 100

# Grab repo directory whether on the cluster or computer
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

results_dir <- paste0(repo_dir, "results/",
                      exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")

aggregated_results <- readRDS(file = paste0(results_dir, "aggregated_results.RDS"))

correlation <- aggregated_results$correlation
metrics <- aggregated_results$metrics
predictions <- aggregated_results$predictions


# Add mean absolute correlation for the final figures for this
mean_abs_cor <- 
  correlation %>% 
  group_by(method, delta, gps_mod, sample_size, sim) %>% 
  dplyr::summarize(pre_cor = mean(pre_cor), post_cor = mean(post_cor)) %>% 
  mutate(covariate = "mean") %>% 
  ungroup()

# Add mean to correlation table
correlation <- 
  rbind(correlation, mean_abs_cor) %>% 
  mutate(covariate = factor(covariate, levels = rev(c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "mean"))))
  
gg_correlation_ipw <-
  correlation %>% 
  filter(method == "ipw") %>% 
  pivot_longer(c("pre_cor", "post_cor"), names_to = "correlation") %>% 
  group_by(gps_mod, covariate, correlation) %>% 
  dplyr::summarize(mean = mean(value),
                   lower = quantile(value, 0.05), 
                   upper = quantile(value, 0.95)) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  ggplot(aes(x = covariate, y = mean, ymin = lower, ymax = upper, color = correlation, group = correlation)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_errorbar(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7)) + 
  #geom_line(size = 0.3) +
  theme_bw() +
  coord_flip() + 
  facet_wrap(~ gps_mod) + 
  labs(y = "Absolute correlation", x = "Covariate", title = "Covariate balance using IPW and entropy based weights")

gg_correlation_default <-
  correlation %>% 
  filter(method == "causal_gps_default") %>% 
  pivot_longer(c("pre_cor", "post_cor"), names_to = "correlation") %>% 
  group_by(gps_mod, covariate, correlation) %>% 
  dplyr::summarize(mean = mean(value),
                   lower = quantile(value, 0.05), 
                   upper = quantile(value, 0.95)) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  ggplot(aes(x = covariate, y = mean, ymin = lower, ymax = upper, color = correlation, group = correlation)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_errorbar(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7)) + 
  #geom_line(size = 0.3) +
  theme_bw() +
  coord_flip() + 
  facet_wrap(~ gps_mod) + 
  labs(y = "Absolute correlation", x = "Covariate", title = "Covariate balance using default causalGPS package")
gg_correlation_default

gg_correlation_tuned <-
  correlation %>% 
  filter(method == "causal_gps_tuned") %>% 
  pivot_longer(c("pre_cor", "post_cor"), names_to = "correlation") %>% 
  group_by(gps_mod, covariate, correlation) %>% 
  dplyr::summarize(mean = mean(value),
                   lower = quantile(value, 0.05), 
                   upper = quantile(value, 0.95)) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  ggplot(aes(x = covariate, y = mean, ymin = lower, ymax = upper, color = correlation, group = correlation)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_errorbar(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7)) + 
  #geom_line(size = 0.3) +
  theme_bw() +
  coord_flip() + 
  facet_wrap(~ gps_mod) +  
  labs(y = "Absolute correlation", x = "Covariate", title = "Covariate balance using tuned causalGPS package")

gg_correlation_tuned

correlation_plot_list <- list(gg_correlation_ipw, gg_correlation_tuned, gg_correlation_default)

# Save ggplot file with all the correlations
#png(filename = paste0(results_dir, "covariate_balance_plots.png"), width = 7, height = 7, units = "in", res = 300)
pdf(file = paste0(results_dir, "covariate_balance_plots.pdf"), width = 8, height = 6)
for (i in 1:length(correlation_plot_list)) {
  print(correlation_plot_list[[i]])
}
dev.off()


# Format data to create plots
plot_data <- 
  metrics %>% 
  mutate(absolute_bias = abs(bias)) %>% 
  tidyr::pivot_longer(c("bias", "absolute_bias", "mse")) %>% 
  mutate(model = factor(model, levels = c("linear_model", "linear_gps", "gam_model", "gam_gps", "causal_gps_default", "causal_gps_tuned", "change_model", "propensity_change")))  %>% 
  dplyr::group_by(model, gps_mod, sample_size, name) %>% 
  dplyr::summarize(mean = mean(value), 
                   lower = quantile(value, 0.1),
                   upper = quantile(value, 0.9)) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction")))

plot_bias <- 
  plot_data %>% 
  filter(model != "causal_gps_default") %>% 
  filter(name == "bias") %>% 
  ggplot() + 
  geom_pointrange(aes(y = mean, ymin = lower, ymax = upper, x = factor(gps_mod), color = model), position = position_dodge(.6)) +
  geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
  labs(x = "Scenario", y = "Value") + 
  theme_bw() + 
  coord_cartesian(ylim = c(-15, 15)) + 
  labs(y = "Bias", x = "Confounding Scenario", title = paste("Bias with", exp_relationship, "relationship"))

plot_abs_bias <- 
  plot_data %>% 
  filter(model != "causal_gps_default") %>%
  filter(name == "absolute_bias") %>% 
  ggplot() + 
  geom_pointrange(aes(y = mean, ymin = lower, ymax = upper, x = factor(gps_mod), color = model), position = position_dodge(.6)) +
  geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
  labs(x = "Scenario", y = "Value") + 
  coord_cartesian(ylim = c(0, 40)) + 
  theme_bw() + 
  labs(y = "Absolute Bias", x = "Confounding Scenario", title = paste("Abolsute bias with", exp_relationship, "relationship"))

plot_mse <- 
  plot_data %>% 
  filter(model != "causal_gps_default") %>% 
  filter(name == "mse") %>% 
  ggplot() + 
  geom_pointrange(aes(y = mean, ymin = lower, ymax = upper, x = factor(gps_mod), color = model), position = position_dodge(.6)) +
  geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype = "dashed", size = 0.2) + 
  labs(x = "Scenario", y = "Value") + 
  theme_bw() + 
  scale_y_continuous(trans='log10') + 
  #coord_cartesian(ylim = c(0, 10)) + 
  labs(y = "Mean squared error (log scale)", x = "Confounding Scenario", title = paste("MSE with", exp_relationship, "relationship"))

bias_plot_list <- list(plot_bias, plot_abs_bias, plot_mse)
# Save ggplot  with bias, absolute bias, mse
#png(filename = paste0(results_dir, "covariate_balance_plots.png"), width = 7, height = 7, units = "in", res = 300)
pdf(file = paste0(results_dir, "bias_mse_plots.pdf"), width = 8, height = 6)
for (i in 1:length(bias_plot_list)) {
  print(bias_plot_list[[i]])
}
dev.off()


# Now make trace plots 
# Make prediction plots
pred_plot <-
  predictions %>%
  dplyr::select(-causal_gps_default) %>% 
  tidyr::pivot_longer(c("linear_model", "linear_gps", "gam_model", "gam_gps", "change_model", "causal_gps_tuned", "true_fit")) %>%
  mutate(name = factor(name, levels = c("linear_model", "linear_gps", "gam_model", "gam_gps", "causal_gps_tuned", "change_model", "propensity_change", "true_fit"))) %>% 
  data.table()

pred_summary <-
  pred_plot[,.(mean = mean(value),
               lower = quantile(value, 0.1),
               upper = quantile(value, 0.9)), by=.(gps_mod, exposure, sample_size, name)] %>% 
  mutate(name = relevel(name, ref = "linear_model")) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  data.table()

# Get a sample of predictions to plot
pred_sample <- 
  pred_plot[sim %in% sample(1:20, 10, replace = F)] %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  data.table()

true_fit_plot <- pred_summary[name == "true_fit"] %>% dplyr::select(-name) %>% data.table()



gg_trace_plot <- 
  ggplot() +
  geom_line(data = pred_summary[name != "true_fit" & sample_size == 1000], aes(x = exposure, y = mean), linetype = "solid") +
  geom_line(data = true_fit_plot[sample_size == 1000], aes(x = exposure, y = mean), color = "black", linetype = "dashed") +
  geom_line(data = pred_sample[name != "true_fit" & sample_size == 1000], aes(x = exposure, y = value, color = factor(sim)), linetype = "solid", alpha = 0.5) + 
  scale_color_manual(values = c("gray", "gray","gray", "gray", "gray", "gray", "gray", "gray", "gray", "gray")) + 
  labs(x = "Exposure", y = "Response") +
  theme_bw() +
  coord_cartesian(ylim = c(-20, 20)) + 
  facet_grid(name ~ gps_mod) + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 16))  

pdf(file = paste0(results_dir, "trace.pdf"), width = 12, height = 10)
print(gg_trace_plot)
dev.off()

# Make an example plot

gg_example <- 
  data.table(exposure = seq(0, 20, length.out = 1000)) %>% 
  mutate(linear = 0.1 * exposure, sublinear = log10(exposure + 1), threshold = ifelse(exposure < 5, 0, 0.1* (exposure - 5))) %>% 
  pivot_longer(-exposure, names_to = "Relationship") %>% 
  ggplot(aes(x = exposure, y = value, color = Relationship)) + 
  geom_line() + 
  labs(x = "Exposure", y = "Outcome") + 
  theme_classic() + 
  theme(text = element_text(size = 18))   

pdf(file = paste0(results_dir, "example.pdf"), width = 7, height = 6)
print(gg_example)
dev.off()
         
