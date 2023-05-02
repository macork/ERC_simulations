# Load packages
rm(list = ls())
library(data.table)
library(tidyverse)
library(latex2exp)
library("CausalGPS", lib.loc = "/n/home_fasse/mcork/apps/ERC_simulation/R_4.0.5")

# Set directory --------
# Grab repo directory whether on the cluster or computer
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "/n/dominici_nsaph_l3/projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

# set figures directory
figure_dir <- paste0(repo_dir, "figures/")

# Load data for each model selection ------------------------------------------------
# First load the linear model example
exp_relationship = "linear"
adjust_confounder = T
time_stamp = "first_draft2"
replicates = 100

results_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")

aggregated_results <- readRDS(file = paste0(results_dir, "aggregated_results.RDS"))

correlation_linear <- aggregated_results$correlation
metrics_linear <- aggregated_results$metrics
predictions_linear <- aggregated_results$predictions

# Change name of models for simplicity 
model_types = c("linear_model", "linear_gps", "gam_model", "gam_gps", "change_model", "change_gps", "causal_gps_default", "causal_gps_tuned")

# collapse by mean bias at each point of exposure
metrics_collapse <- 
  metrics_linear[, .(abs_bias = abs(mean(bias)),
                     mse = mean(mse),
                     rmse = sqrt(mean(mse))), 
                 by = .(model, exposure, gps_mod, sample_size, confounder_setting, out_relationship)]

metrics_collapse_linear <- 
  metrics_collapse[, .(abs_bias = mean(abs_bias),
                       mse = mean(mse),
                       rmse = mean(rmse)), by = .(model, gps_mod, sample_size, confounder_setting, out_relationship)] %>% 
  mutate(model = factor(model, levels = model_types))  %>% 
  mutate(model  = fct_recode(model, linear = "linear_model", linear_ent = "linear_gps", gam = "gam_model", 
                             gam_ent = "gam_gps", change = "change_model", change_ent = "change_gps",
                             causalGPS = "causal_gps_tuned")) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  mutate(out_relationship = factor(out_relationship, levels = c("linear", "interaction"))) %>% 
  mutate(confounder_setting = factor(confounder_setting, levels = c("simple", "nonzero")))

# Add correct labels  
levels(metrics_collapse_linear$out_relationship) = c(linear = TeX("$\\mu_{linear}$"), interaction = TeX("$\\mu_{linear, int}$"))


# Rename predictions to match with new model names
predictions_linear <- 
  predictions_linear %>% 
  rename(linear = linear_model, linear_ent = linear_gps, gam = gam_model,
         gam_ent = gam_gps, causalGPS = causal_gps_tuned, change = change_model, change_ent = change_gps) %>% 
  data.table()


# load sublinear ----------------------------------------------------------------------------------------
exp_relationship = "sublinear"
adjust_confounder = T
time_stamp = "first_draft"
replicates = 100

results_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")

aggregated_results <- readRDS(file = paste0(results_dir, "aggregated_results.RDS"))

correlation_sublinear <- aggregated_results$correlation
metrics_sublinear <- aggregated_results$metrics
predictions_sublinear <- aggregated_results$predictions

# collapse by mean bias at each point of exposure
metrics_collapse <- 
  metrics_sublinear[, .(abs_bias = abs(mean(bias)),
                        mse = mean(mse),
                        rmse = sqrt(mean(mse))), 
                    by = .(model, exposure, gps_mod, sample_size, confounder_setting, out_relationship)]

metrics_collapse_sublinear <- 
  metrics_collapse[, .(abs_bias = mean(abs_bias),
                       mse = mean(mse),
                       rmse = mean(rmse)), by = .(model, gps_mod, sample_size, confounder_setting, out_relationship)] %>% 
  mutate(model = factor(model, levels = model_types))  %>% 
  mutate(model  = fct_recode(model, linear = "linear_model", linear_ent = "linear_gps", gam = "gam_model", 
                             gam_ent = "gam_gps", change = "change_model", change_ent = "change_gps",
                             causalGPS = "causal_gps_tuned")) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  mutate(out_relationship = factor(out_relationship, levels = c("linear", "interaction"))) %>% 
  mutate(confounder_setting = factor(confounder_setting, levels = c("simple", "nonzero")))

# Add correct labels  
levels(metrics_collapse_sublinear$out_relationship) = c(linear = TeX("$\\mu_{sublinear}$"), interaction = TeX("$\\mu_{sublinear, int}$"))

predictions_sublinear <- 
  predictions_sublinear %>%
  rename(linear = linear_model, linear_ent = linear_gps, gam = gam_model,
         gam_ent = gam_gps, causalGPS = causal_gps_tuned, change = change_model, change_ent = change_gps) %>% 
  data.table()

#load threshold ----------------------------------------------------------------------------------------
exp_relationship = "threshold"
adjust_confounder = T
time_stamp = "first_draft2"
replicates = 100

results_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")

aggregated_results <- readRDS(file = paste0(results_dir, "aggregated_results.RDS"))

correlation_threshold <- aggregated_results$correlation
metrics_threshold <- aggregated_results$metrics
predictions_threshold <- aggregated_results$predictions

# collapse by mean bias at each point of exposure
metrics_collapse <- 
  metrics_threshold[, .(abs_bias = abs(mean(bias)),
                        mse = mean(mse),
                        rmse = sqrt(mean(mse))), 
                    by = .(model, exposure, gps_mod, sample_size, confounder_setting, out_relationship)]

metrics_collapse_threshold <- 
  metrics_collapse[, .(abs_bias = mean(abs_bias),
                       mse = mean(mse),
                       rmse = mean(rmse)), by = .(model, gps_mod, sample_size, confounder_setting, out_relationship)] %>% 
  mutate(model = factor(model, levels = model_types))  %>% 
  mutate(model  = fct_recode(model, linear = "linear_model", linear_ent = "linear_gps", gam = "gam_model", 
                             gam_ent = "gam_gps", change = "change_model", change_ent = "change_gps",
                             causalGPS = "causal_gps_tuned")) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  mutate(out_relationship = factor(out_relationship, levels = c("linear", "interaction"))) %>% 
  mutate(confounder_setting = factor(confounder_setting, levels = c("simple", "nonzero")))

# Add correct labels  
levels(metrics_collapse_threshold$out_relationship) = c(linear = TeX("$\\mu_{threshold}$"), interaction = TeX("$\\mu_{threshold, int}$"))


predictions_threshold <- 
  predictions_threshold %>% 
  rename(linear = linear_model, linear_ent = linear_gps, gam = gam_model,
         gam_ent = gam_gps, causalGPS = causal_gps_tuned, change = change_model, change_ent = change_gps) %>% 
  data.table()


# Now make figures 1-3

# Calculate mean absolute correlation
mean_abs_cor <- 
  correlation_linear %>% 
  pivot_longer(c("pre_cor", "post_cor"), names_to = "correlation") %>% 
  group_by(gps_mod, correlation, sample_size, confounder_setting, out_relationship, method, sim) %>% 
  dplyr::summarize(value = mean(value)) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  mutate(covariate = "mean") %>% 
  data.table()

balance_figure_1 <- 
  mean_abs_cor %>% 
  filter(correlation == "post_cor", out_relationship == "interaction", confounder_setting == "simple") %>%
  filter(method %in% c("causal_gps_tuned", "ipw")) %>% 
  mutate(method = gsub("causal_gps_tuned", "CausalGPS", method)) %>% 
  mutate(method = gsub("ipw", "Entropy", method)) %>% 
  group_by(gps_mod, correlation, sample_size, confounder_setting, out_relationship, method) %>% 
  dplyr::summarize(balance_pct = mean(value < 0.1)) %>% 
  ungroup() %>% 
  mutate(confounder_setting = factor(confounder_setting, levels = c("simple", "nonzero"))) %>% 
  data.table() %>% 
  ggplot(aes(x = gps_mod, y = balance_pct, color = method, group = method)) +
  geom_point() + 
  geom_line(size = 0.3) +
  theme_bw(base_size = 14) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size = 12)) + 
  theme(
    legend.position = c(0.85, 0.2),
    legend.title = element_blank(),
    legend.text=element_text(size=16)
  ) + 
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  facet_wrap(~ sample_size) + 
  scale_color_discrete("") + 
  labs(y = "Proportion of simulations that\nmeet balance threshold", x = "Exposure model") # 

ggsave(balance_figure_1, file = paste0(figure_dir, "fig1_balance.png"),
       width = 10, height = 4,  dpi = 400)

# Put together data for correlation figure
correlation_linear_summary <- 
  correlation_linear %>% 
  # filter(sample_size == 1000, confounder_setting == "simple", out_relationship == "linear", method == "ipw") %>% 
  pivot_longer(c("pre_cor", "post_cor"), names_to = "correlation") %>% 
  dplyr::select(gps_mod, covariate, correlation, sample_size, confounder_setting, out_relationship, method, sim, value) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  rbind(mean_abs_cor) %>% 
  group_by(gps_mod, covariate, correlation, sample_size, confounder_setting, out_relationship, method) %>% 
  dplyr::summarize(mean = mean(value),
                   lower = quantile(value, 0.025), 
                   upper = quantile(value, 0.975)) %>%
  mutate(covariate = factor(covariate, levels = rev(c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "mean")))) %>% 
  data.table()

# Look at simple confounder setting IPW
correlation_plot_entropy <- 
  correlation_linear_summary %>% 
  filter(sample_size == 1000, confounder_setting == "simple", out_relationship == "linear", method == "ipw") %>% 
  ggplot(aes(x = covariate, y = mean, ymin = lower, ymax = upper, color = correlation, group = correlation)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_errorbar(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7)) + 
  scale_color_discrete("", breaks = c("pre_cor", "post_cor"), labels=c("Pre-weighting", "Post-weighting")) + 
  #geom_line(size = 0.3) +
  theme_bw(base_size = 14) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size = 12)) + 
  theme(legend.position = "bottom") + 
  coord_flip() + 
  facet_wrap(~ gps_mod) + 
  labs(y = "Absolute correlation", x = "Covariate")

ggsave(correlation_plot_entropy, file = paste0(figure_dir, "fig2_correlation_plot_entropy.png"),
       width = 10, height = 7,  dpi = 400)

# Repeat for CausalGPS package
correlation_plot_CausalGPS <- 
  correlation_linear_summary %>% 
  filter(sample_size == 1000, confounder_setting == "simple", out_relationship == "linear", method == "causal_gps_tuned") %>% 
  ggplot(aes(x = covariate, y = mean, ymin = lower, ymax = upper, color = correlation, group = correlation)) +
  geom_hline(aes(yintercept = 0.1)) + 
  geom_errorbar(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7)) + 
  scale_color_discrete("", breaks = c("pre_cor", "post_cor"), labels=c("Pre-weighting", "Post-weighting")) + 
  theme_bw(base_size = 14) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size = 12)) + 
  theme(legend.position = "bottom") + 
  coord_flip() + 
  facet_wrap(~ gps_mod) + 
  labs(y = "Absolute correlation", x = "Covariate")

ggsave(correlation_plot_CausalGPS, file = paste0(figure_dir, "fig3_correlation_plot_causalGPS.png"),
       width = 10, height = 7,  dpi = 400)


# Now work on bias and RMSE plots 
# Function to make bias and RMSE plots for outcome model setting
bias_RMSE_plots <- function(out_model) {
  
  if (out_model == "linear") {
    metrics_collapse = metrics_collapse_linear
  } else if (out_model == "sublinear") {
    metrics_collapse = metrics_collapse_sublinear
  } else if (out_model == "threshold") {
    metrics_collapse = metrics_collapse_threshold
  } else {
    stop("Outcome model must be linear, sublinear or threshold")
  }
  
  # Absolute bias plot
  abs_bias_plot <- 
    metrics_collapse %>% 
    filter(model != "causal_gps_default") %>% 
    filter(confounder_setting == "simple") %>% 
    #filter(sample_size == 1000) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
    geom_point(aes(y = abs_bias, x = factor(gps_mod), 
                   color = model, fill = model, shape = model), position = position_dodge(.85), size = 1.8) + 
    geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
    labs(x = "Scenario", y = "Value") + 
    theme_bw(base_size = 14) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 1)) + 
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 12)) + 
    scale_shape_manual(values=c(21, 10, 22, 12, 23, 9, 8))+
    #coord_cartesian(ylim = c(-15, 15)) + 
    facet_grid(sample_size ~ out_relationship, scales = "free", labeller = label_parsed) + 
    #facet_wrap(~ out_relationship, scales = "free") + 
    labs(y = "Absolute bias", x = "Exposure model")
  
  # Now plot the RMSE 
  rmse_plot <- 
    metrics_collapse %>% 
    filter(model != "causal_gps_default") %>% 
    filter(confounder_setting == "simple") %>% 
    #filter(sample_size == 1000) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
    geom_point(aes(y = rmse, x = factor(gps_mod), color = model, fill = model, shape = model), position = position_dodge(.85), size = 1.8) + 
    # geom_pointrange(aes(y = mean, ymin = mean - sd, ymax = mean + sd, x = factor(gps_mod), color = model), position = position_dodge(.85), size = 0.35) +
    geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
    labs(x = "Scenario", y = "Value") + 
    # scale_y_continuous(trans='log10') + 
    theme_bw(base_size = 14) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 1)) + 
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 12)) + 
    scale_shape_manual(values=c(21, 10, 22, 12, 23, 9, 8)) +
    #coord_cartesian(ylim = c(-15, 15)) + 
    facet_grid(sample_size ~ out_relationship, scales = "free", labeller = label_parsed) + 
    #facet_wrap(~ out_relationship, scales = "free") + 
    labs(y = "Root mean squared error (RMSE)", x = "Exposure model")
  
  # Now return both plots
  return(list(bias_plot = abs_bias_plot, rmse_plot = rmse_plot))
}

linear_plots <- bias_RMSE_plots("linear")
abs_bias_linear <- linear_plots[[1]]
rmse_linear <- linear_plots[[2]]

# Save plots
ggsave(
  filename = paste0(figure_dir, "fig4_abs_bias_linear.png"),
  plot = abs_bias_linear,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  filename = paste0(figure_dir, "fig5_rmse_linear.png"),
  plot = rmse_linear,
  width = 10,
  height = 7,
  dpi = 400
)

# Now for sublinear ----------------------------------------
sublinear_plots <- bias_RMSE_plots("sublinear")
abs_bias_sublinear <- sublinear_plots[[1]]
rmse_sublinear <- sublinear_plots[[2]]

# Save plots
ggsave(
  filename = paste0(figure_dir, "fig7_abs_bias_sublinear.png"),
  plot = abs_bias_sublinear,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  filename = paste0(figure_dir, "fig8_rmse_sublinear.png"),
  plot = rmse_sublinear,
  width = 10,
  height = 7,
  dpi = 400
)


# Now for threshold ----------------------------------------
threshold_plots <- bias_RMSE_plots("threshold")
abs_bias_threshold<- threshold_plots[[1]]
rmse_threshold <- threshold_plots[[2]]

# Save plots
ggsave(
  filename = paste0(figure_dir, "fig10_abs_bias_threshold.png"),
  plot = abs_bias_threshold,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  filename = paste0(figure_dir, "fig11_rmse_threshold.png"),
  plot = rmse_threshold,
  width = 10,
  height = 7,
  dpi = 400
)

# Now work on plots of ERC under heterogeneous treatment setting ----------------------------------------
ERC_plot <- function(out_model) {
  
  if (out_model == "linear") {
    preds = predictions_linear
  } else if (out_model == "sublinear") {
    preds = predictions_sublinear
  } else if (out_model == "threshold") {
    preds = predictions_threshold
  } else {
    stop("Outcome model must be linear, sublinear or threshold")
  }
  
  # Prediction plots for ERC under heterogeneous treatment effects
  pred_plot <-
    preds %>%
    filter(confounder_setting == "simple", out_relationship == "interaction", sample_size == 1000) %>% 
    dplyr::select(-causal_gps_default) %>% 
    tidyr::pivot_longer(c("linear", "linear_ent", "gam", "gam_ent", "change","change_ent", "causalGPS", "true_fit")) %>%
    mutate(name = factor(name, levels = c("linear", "linear_ent", "gam", "gam_ent", "change", "change_ent", "causalGPS", "true_fit"))) %>% 
    data.table()
  
  # Now create summary statistic based on mean/lower/upper 
  pred_summary <-
    pred_plot[,.(mean = mean(value),
                 lower = quantile(value, 0.05),
                 upper = quantile(value, 0.95)), by=.(exposure, gps_mod, sample_size, confounder_setting, out_relationship, name)] %>% 
    mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                      gps_mod == 2 ~ "heavytail", 
                                      gps_mod == 3 ~ "nonlinear", 
                                      gps_mod == 4 ~ "interaction"), levels = c("linear", "heavytail", "nonlinear", "interaction"))) %>% 
    data.table()
  
  # Get a sample of predictions to plot
  pred_sample <- 
    pred_plot[sim %in% sample(1:20, 10, replace = F)] %>% 
    mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                      gps_mod == 2 ~ "heavytail", 
                                      gps_mod == 3 ~ "nonlinear", 
                                      gps_mod == 4 ~ "interaction"), levels = c("linear", "heavytail", "nonlinear", "interaction"))) %>% 
    data.table()
  
  levels(pred_summary$gps_mod) <- c(linear = TeX("$E_{linear}$"), heavytail = TeX("$E_{heavytail}$"), nonlinear = TeX("$E_{nonlinear}$"), interaction = TeX("$E_{interaction}$"))
  levels(pred_sample$gps_mod) <- c(linear = TeX("$E_{linear}$"), heavytail = TeX("$E_{heavytail}$"), nonlinear = TeX("$E_{nonlinear}$"), interaction = TeX("$E_{interaction}$"))
  
  # Plot the true fit
  true_fit_plot <- pred_summary[name == "true_fit"] %>% dplyr::select(-name) %>% data.table()
  
  gg_ERC <- 
    ggplot() +
    geom_line(data = pred_sample[name != "true_fit" & sample_size == 1000], aes(x = exposure, y = value, color = factor(sim)), linetype = "solid", alpha = 0.5, size = 0.3) + 
    scale_color_manual(values = c("gray", "gray","gray", "gray", "gray", "gray", "gray", "gray", "gray", "gray")) + 
    geom_line(data = pred_summary[name != "true_fit" & sample_size == 1000], aes(x = exposure, y = mean), linetype = "solid", size = 0.3) +
    geom_line(data = true_fit_plot[sample_size == 1000], aes(x = exposure, y = mean), color = "black", linetype = "dashed", size = 0.3) +
    labs(x = "Exposure", y = "Response", title = "") +
    theme_bw(base_size = 14) + 
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 12)) + 
    coord_cartesian(ylim = c(10, 50)) + 
    facet_grid(name ~ gps_mod, labeller = label_parsed) + 
    theme(legend.position = "none")
  return(gg_ERC)
}

linear_ERC <- ERC_plot("linear")

# Save plots to use in results
ggsave(
  filename = paste0(figure_dir, "fig6_ERC_linear.png"),
  plot = linear_ERC,
  width = 7,
  height = 10,
  dpi = 400
)

sublinear_ERC <- ERC_plot("sublinear")

# Save plots to use in results
ggsave(
  filename = paste0(figure_dir, "fig9_ERC_sublinear.png"),
  plot = sublinear_ERC,
  width = 7,
  height = 10,
  dpi = 400
)

threshold_ERC <- ERC_plot("threshold")

# Save plots to use in results
ggsave(
  filename = paste0(figure_dir, "fig12_ERC_threshold.png"),
  plot = threshold_ERC,
  width = 7,
  height = 10,
  dpi = 400
)


#### Last 3 figures are produced in data_application/2_post_processing.R code 




