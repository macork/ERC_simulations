# Load packages
rm(list = ls())
library(tidyverse)
library(latex2exp)
library(CausalGPS)
library(patchwork)


#unlink("/n/home_fasse/mcork/apps/ERC_simulation/R_4.0.5/00LOCK-vctrs", recursive = T)
#install.packages("tidyverse", dependencies = T)

# Set directory --------
# Grab repo directory (based on either local or on cluster)
if (Sys.getenv("USER") == "mcork") {
  repo_dir <- "~/nsaph_projects/ERC_simulation/Simulation_studies/"
} else if (Sys.getenv("USER") == "michaelcork") {
  repo_dir <- "~/Desktop/Francesca_research/Simulation_studies/"
}

# set figures directory
figure_dir <- paste0(repo_dir, "figures_resubmission/")

# Function to summarize metrics
collapse_metrics <- function(metrics_type, model_types) {
  
  # collapse by mean bias at each point of exposure
  metrics_collapse <- 
    metrics_type %>%
    filter(model %in% model_types) %>%
    group_by(model, exposure, gps_mod, sample_size, outcome_interaction) %>%
    summarize(
      abs_bias = abs(mean(bias, na.rm = TRUE)),
      mse = mean(mse, na.rm = TRUE),
      rmse = sqrt(mean(mse, na.rm = TRUE)),
      .groups = 'drop'
    )
  
  metrics_collapse <- 
    metrics_collapse %>%
    group_by(model, gps_mod, sample_size, outcome_interaction) %>%
    summarize(
      abs_bias = mean(abs_bias, na.rm = TRUE),
      mse = mean(mse, na.rm = TRUE),
      rmse = mean(rmse, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Model = factor(model, 
                     levels = c("linear_model", "ent_linear", "gam_model", "ent_gam", 
                                "change_model", "ent_change", "causal_gps_tuned"), 
                     labels = c("Linear", "Linear\nentropy", "GAM", "GAM\nentropy", 
                                "Change", "Change\nentropy", "CausalGPS")), 
      gps_mod = factor(case_when(
        gps_mod == 1 ~ "linear", 
        gps_mod == 2 ~ "heavy tail", 
        gps_mod == 3 ~ "nonlinear", 
        gps_mod == 4 ~ "interaction"), 
        levels = c("linear", "heavy tail", "nonlinear", "interaction")), 
      out_relationship = if_else(as.logical(outcome_interaction), "interaction", "linear"), 
      out_relationship = factor(out_relationship, levels = c("linear", "interaction"))
    )
  
  # Add correct labels  
  levels(metrics_collapse$out_relationship) = c(linear = TeX("$\\mu_{linear}$"), interaction = TeX("$\\mu_{linear, int}$"))
  
  return(metrics_collapse)
}


# Load data for each model selection ------------------------------------------------
# First load the linear model example
exp_relationship = "linear"
adjust_confounder = T
#time_stamp = "first_draft2"
#time_stamp = "first_submission4" # With trimming in CausalGPS
time_stamp = "resubmission_thousand"
replicates = 1000

results_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")

aggregated_results <- readRDS(file = paste0(results_dir, "aggregated_results.RDS"))

correlation_linear <- aggregated_results$correlation
metrics_linear <- aggregated_results$metrics
predictions_linear <- aggregated_results$predictions
convergence_linear <- aggregated_results$convergence

# Change name of models for simplicity 
model_types = c("linear_model", "ent_linear",
                "gam_model", "ent_gam", 
                "change_model", "ent_change",
                "causal_gps_tuned")

# collapse by mean bias at each point of exposure
metrics_collapse_linear <- collapse_metrics(metrics_linear, model_types)

# Rename predictions to match with new model names
predictions_linear <- 
  predictions_linear %>% 
  rename(linear = linear_model, linear_ent = ent_linear, gam = gam_model,
         gam_ent = ent_gam, CausalGPS = causal_gps_tuned, change = change_model, change_ent = ent_change) %>% 
  mutate(out_relationship = if_else(as.logical(outcome_interaction), "interaction", "linear"), 
         out_relationship = factor(out_relationship, levels = c("linear", "interaction")))


# load sublinear ----------------------------------------------------------------------------------------
exp_relationship = "sublinear"
adjust_confounder = T
#time_stamp = "first_submission4"
time_stamp = "resubmission_thousand"
replicates = 1000

results_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")

aggregated_results <- readRDS(file = paste0(results_dir, "aggregated_results.RDS"))

correlation_sublinear <- aggregated_results$correlation
metrics_sublinear <- aggregated_results$metrics
predictions_sublinear <- aggregated_results$predictions
convergence_sublinear <- aggregated_results$convergence

# collapse by mean bias at each point of exposure
metrics_collapse_sublinear <- collapse_metrics(metrics_sublinear, model_types)

# Rename predictions to match with new model names
predictions_sublinear <- 
  predictions_sublinear %>% 
  rename(linear = linear_model, linear_ent = ent_linear, gam = gam_model,
         gam_ent = ent_gam, CausalGPS = causal_gps_tuned, change = change_model, change_ent = ent_change) %>% 
  mutate(out_relationship = if_else(as.logical(outcome_interaction), "interaction", "linear"), 
         out_relationship = factor(out_relationship, levels = c("linear", "interaction")))


#load threshold ----------------------------------------------------------------------------------------
exp_relationship = "threshold"
adjust_confounder = T
#time_stamp = "first_submission4"
time_stamp = "resubmission_thousand"
replicates = 1000

results_dir <- paste0(repo_dir, "results/", exp_relationship, "_", adjust_confounder, "/", time_stamp, "/")

aggregated_results <- readRDS(file = paste0(results_dir, "aggregated_results.RDS"))

correlation_threshold <- aggregated_results$correlation
metrics_threshold <- aggregated_results$metrics
predictions_threshold <- aggregated_results$predictions

# collapse by mean bias at each point of exposure
metrics_collapse_threshold <- collapse_metrics(metrics_threshold, model_types)

predictions_threshold <- 
  predictions_threshold %>% 
  rename(linear = linear_model, linear_ent = ent_linear, gam = gam_model,
         gam_ent = ent_gam, CausalGPS = causal_gps_tuned, change = change_model, change_ent = ent_change) %>% 
  mutate(out_relationship = if_else(as.logical(outcome_interaction), "interaction", "linear"), 
         out_relationship = factor(out_relationship, levels = c("linear", "interaction")))

# Now make figures 1 and 2 for supplementary materials ------------------------------------------------------

# Calculate mean absolute correlation of all six covariates
mean_abs_cor <- 
  correlation_linear %>% 
  filter(method %in% c("ent", "causal_gps_tuned")) %>% 
  pivot_longer(cols = c("pre_cor", "post_cor"), names_to = "correlation") %>% 
  group_by(gps_mod, correlation, sample_size, outcome_interaction, method, sim) %>% 
  summarize(value = mean(value, na.rm = TRUE), .groups = 'drop') %>% 
  mutate(
    gps_mod = factor(gps_mod, levels = 1:4, labels = c("linear", "heavy tail", "nonlinear", "interaction")),
    covariate = "mean"
  )

balance_figure_1 <- 
  mean_abs_cor %>% 
  filter(correlation == "post_cor", outcome_interaction == "T") %>%
  filter(method %in% c("causal_gps_tuned", "ent")) %>% 
  mutate(
    method = case_when(
      method == "causal_gps_tuned" ~ "CausalGPS",
      method == "ent" ~ "Entropy",
      TRUE ~ method)
    ) %>%
  group_by(gps_mod, correlation, sample_size, outcome_interaction, method) %>% 
  summarize(balance_pct = mean(value < 0.1), .groups = 'drop') %>%
  ggplot(aes(x = gps_mod, y = balance_pct, color = method, group = method)) +
  geom_point() + 
  geom_line(size = 0.3) +
  theme_bw(base_size = 14) + 
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    legend.position = c(0.85, 0.2),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  ) + 
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  facet_wrap(~ sample_size) + 
  scale_color_discrete("") + 
  labs(y = "Proportion of simulations that\nmeet balance threshold", 
       x = "Exposure model")


ggsave(balance_figure_1, file = paste0(figure_dir, "supp_fig1_balance.png"),
       width = 10, height = 4,  dpi = 400)

# Put together data for correlation figure
correlation_linear_summary <- 
  correlation_linear %>% 
  # filter(sample_size == 1000, confounder_setting == "simple", out_relationship == "linear", method == "ipw") %>% 
  pivot_longer(c("pre_cor", "post_cor"), names_to = "correlation") %>% 
  dplyr::select(gps_mod, covariate, correlation, sample_size, outcome_interaction, method, sim, value) %>% 
  mutate(gps_mod = factor(case_when(gps_mod == 1 ~ "linear", 
                                    gps_mod == 2 ~ "heavy tail", 
                                    gps_mod == 3 ~ "nonlinear", 
                                    gps_mod == 4 ~ "interaction"), levels = c("linear", "heavy tail", "nonlinear", "interaction"))) %>% 
  rbind(mean_abs_cor) %>% 
  group_by(gps_mod, covariate, correlation, sample_size, outcome_interaction, method) %>% 
  dplyr::summarize(mean = mean(value),
                   lower = quantile(value, 0.025), 
                   upper = quantile(value, 0.975)) %>%
  mutate(covariate = factor(covariate, levels = rev(c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6", "mean")))) 

# Look at simple confounder setting IPW
correlation_plot_entropy <- 
  correlation_linear_summary %>% 
  filter(sample_size == 1000, outcome_interaction == "F", method == "ent") %>% 
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

# ggsave(correlation_plot_entropy, file = paste0(figure_dir, "supp_fig2_correlation_plot_entropy.png"),
#        width = 10, height = 7,  dpi = 400)

# Repeat for CausalGPS package
correlation_plot_CausalGPS <- 
  correlation_linear_summary %>% 
  filter(sample_size == 1000, outcome_interaction == "T", method == "causal_gps_tuned") %>% 
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

# ggsave(correlation_plot_CausalGPS, file = paste0(figure_dir, "supp_fig3_correlation_plot_causalGPS.png"),
#        width = 10, height = 7,  dpi = 400)

# wrap the two plots with brackets
plot_a <- correlation_plot_entropy
plot_b <- correlation_plot_CausalGPS

# arrange the two plots side by side
combined_plots <- plot_a + plot_b & theme(legend.position = "bottom")
combined_plots <- combined_plots + plot_layout(guides = "collect", ncol = 1)

# arrange the two plots side by side
combined_plots <- plot_a + plot_b & theme(legend.position = "bottom")
combined_plots <- 
  combined_plots + 
  plot_layout(guides = "collect", ncol = 1) + 
  plot_annotation(tag_levels = 'a')

ggsave(combined_plots, file = paste0(figure_dir, "supp_fig2_correlation_plot.png"),
       width = 10, height = 14,  dpi = 400)


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
    # filter(model != "causal_default") %>% 
    # filter(confounder_setting == "simple") %>% 
    #filter(sample_size == 1000) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
    geom_point(aes(y = abs_bias, x = factor(gps_mod), 
                   color = Model, fill = Model, shape = Model), position = position_dodge(.85), size = 1.8) + 
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
    # filter(model != "causal_gps_default") %>% 
    # filter(confounder_setting == "simple") %>% 
    #filter(sample_size == 1000) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
    geom_point(aes(y = rmse, x = factor(gps_mod), color = Model, fill = Model, shape = Model), position = position_dodge(.85), size = 1.8) + 
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
  filename = paste0(figure_dir, "supp_fig3_abs_bias_linear.pdf"),
  plot = abs_bias_linear,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  filename = paste0(figure_dir, "supp_fig4_rmse_linear.pdf"),
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
  filename = paste0(figure_dir, "supp_fig6_abs_bias_sublinear.pdf"),
  plot = abs_bias_sublinear,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  filename = paste0(figure_dir, "supp_fig7_rmse_sublinear.pdf"),
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
  filename = paste0(figure_dir, "supp_fig9_abs_bias_threshold.pdf"),
  plot = abs_bias_threshold,
  width = 10,
  height = 7,
  dpi = 400
)

ggsave(
  filename = paste0(figure_dir, "supp_fig10_rmse_threshold.pdf"),
  plot = rmse_threshold,
  width = 10,
  height = 7,
  dpi = 400
)

# Now create plots for figure
bias_RMSE_combine_plot <- function(metrics_collapse) {
  
  # Absolute bias plot
  abs_bias_plot <- 
    metrics_collapse %>% 
    filter(sample_size %in% c(1000)) %>%
    # filter(model != "causal_gps_default") %>% 
    # filter(confounder_setting == "simple") %>% 
    #filter(sample_size == 1000) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
    geom_point(aes(y = abs_bias, x = factor(gps_mod), 
                   color = Model, fill = Model, shape = Model), 
               position = position_dodge(.85), size = 2.7) + 
    geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
    labs(x = "Scenario", y = "Value") + 
    theme_bw(base_size = 14) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 1)) + 
    theme(text = element_text(size=14),
          axis.text.x = element_text(size = 12),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank()) + 
    scale_shape_manual(values=c(21, 10, 22, 12, 23, 9, 8))+
    #coord_cartesian(ylim = c(-15, 15)) + 
    facet_grid( ~ out_relationship, scales = "free", labeller = label_parsed) + 
    #facet_wrap(~ out_relationship, scales = "free") + 
    labs(y = "Absolute bias", x = "")
  
  # Now plot the RMSE 
  rmse_plot <- 
    metrics_collapse %>% 
    filter(sample_size %in% c(1000)) %>%
    # filter(model != "causal_gps_default") %>% 
    # filter(confounder_setting == "simple") %>% 
    #filter(sample_size == 1000) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
    geom_point(aes(y = rmse, x = factor(gps_mod), 
                   color = Model, fill = Model, shape = Model),
               position = position_dodge(.85), size = 2.7) + 
    # geom_pointrange(aes(y = mean, ymin = mean - sd, ymax = mean + sd, x = factor(gps_mod), color = model), position = position_dodge(.85), size = 0.35) +
    geom_vline(xintercept=c(1.5, 2.5,3.5, 4.5, 5.5), linetype="dashed", size = 0.2) + 
    labs(x = "Scenario", y = "Value") + 
    # scale_y_continuous(trans='log10') + 
    theme_bw(base_size = 14) + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 1)) + 
    theme(text = element_text(size=14),
          axis.title.x = element_text(vjust = -1),
          axis.text.x = element_text(size = 12),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank()) + 
    scale_shape_manual(values=c(21, 10, 22, 12, 23, 9, 8)) +
    #coord_cartesian(ylim = c(-15, 15)) + 
    facet_grid( ~ out_relationship, scales = "free", labeller = label_parsed) + 
    #facet_wrap(~ out_relationship, scales = "free") + 
    labs(y = "RMSE", x = "Exposure model")
  
  # wrap the two plots with brackets
  plot_a <- abs_bias_plot
  plot_b <- rmse_plot
  
  # arrange the two plots side by side
  combined_plots <- plot_a + plot_b & theme(legend.position = "bottom")
  combined_plots <- 
    combined_plots + 
    plot_layout(guides = "collect", ncol = 1) + 
    plot_annotation(tag_levels = 'a')
  
  # Now return both plots
  return(combined_plots)
}

linear_combine_plot <- bias_RMSE_combine_plot(metrics_collapse_linear)

# Save plots
ggsave(
  filename = paste0(figure_dir, "fig1_linear_combine_plot.pdf"),
  plot = linear_combine_plot,
  width = 10,
  height = 9,
  dpi = 400
)

sublinear_combine_plot <- bias_RMSE_combine_plot(metrics_collapse_sublinear)

# Save plots
ggsave(
  filename = paste0(figure_dir, "fig2_sublinear_combine_plot.pdf"),
  plot = sublinear_combine_plot,
  width = 10,
  height = 9,
  dpi = 400
)

threshold_combine_plot <- bias_RMSE_combine_plot(metrics_collapse_threshold)

# Save plots
ggsave(
  filename = paste0(figure_dir, "fig3_threshold_combine_plot.pdf"),
  plot = threshold_combine_plot,
  width = 10,
  height = 9,
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
  
  pred_plot <- 
    preds %>%
    filter(out_relationship == "interaction", sample_size == 1000) %>% 
    pivot_longer(cols = c("linear", "linear_ent", "gam", "gam_ent", "change", "change_ent", "CausalGPS", "true_fit"), 
                 names_to = "name") %>%
    mutate(name = factor(name, levels = c("linear", "linear_ent", "gam", "gam_ent", "change", "change_ent", "CausalGPS", "true_fit"))) 
  
  pred_summary <- 
    pred_plot %>%
    group_by(exposure, gps_mod, sample_size, out_relationship, name) %>%
    summarize(
      mean = mean(value, na.rm = TRUE),
      lower = quantile(value, 0.05, na.rm = TRUE),
      upper = quantile(value, 0.95, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(gps_mod = recode_factor(gps_mod, `1` = "linear", `2` = "heavytail", `3` = "nonlinear", `4` = "interaction"))
  
  pred_sample <- 
    pred_plot %>%
    filter(sim %in% sample(1:900, 10, replace = FALSE)) %>%
    mutate(gps_mod = recode_factor(gps_mod, `1` = "linear", `2` = "heavytail", `3` = "nonlinear", `4` = "interaction"))
  
  adjust_gps_levels <- function(data) {
    levels(data$gps_mod) <- c(linear = TeX("$E_{linear}$"), heavytail = TeX("$E_{heavytail}$"), nonlinear = TeX("$E_{nonlinear}$"), interaction = TeX("$E_{interaction}$"))
    data
  }
  
  pred_summary <- adjust_gps_levels(pred_summary)
  pred_sample <- adjust_gps_levels(pred_sample)
  
  true_fit_plot <- 
    pred_summary %>%
    filter(name == "true_fit") %>%
    select(-name)
  
  gg_ERC <- 
    ggplot() +
    geom_line(data = filter(pred_sample, name != "true_fit", sample_size == 1000), 
              aes(x = exposure, y = value, color = factor(sim)), 
              linetype = "solid", alpha = 0.5, size = 0.3) + 
    scale_color_manual(values = rep("gray", 10)) + 
    geom_line(data = filter(pred_summary, name != "true_fit", sample_size == 1000), 
              aes(x = exposure, y = mean), 
              linetype = "solid", size = 0.3) +
    geom_line(data = filter(true_fit_plot, sample_size == 1000), 
              aes(x = exposure, y = mean), 
              color = "black", linetype = "dashed", size = 0.3) +
    labs(x = "Exposure", y = "Response", title = "") +
    theme_bw(base_size = 14) + 
    theme(text = element_text(size = 14),
          axis.text.x = element_text(size = 12)) + 
    coord_cartesian(ylim = c(10, 50)) + 
    facet_grid(name ~ gps_mod, labeller = label_parsed) + 
    theme(legend.position = "none")
  
  return(gg_ERC)
}


linear_ERC <- ERC_plot("linear")

# Save plots to use in results
ggsave(
  filename = paste0(figure_dir, "supp_fig5_ERC_linear.png"),
  plot = linear_ERC,
  width = 7,
  height = 10,
  dpi = 400
)

sublinear_ERC <- ERC_plot("sublinear")

# Save plots to use in results
ggsave(
  filename = paste0(figure_dir, "supp_fig8_ERC_sublinear.png"),
  plot = sublinear_ERC,
  width = 7,
  height = 10,
  dpi = 400
)

threshold_ERC <- ERC_plot("threshold")

# Save plots to use in results
ggsave(
  filename = paste0(figure_dir, "supp_fig11_ERC_threshold.png"),
  plot = threshold_ERC,
  width = 7,
  height = 10,
  dpi = 400
)


## Now working on balancing on second moment and how that affects the simulation
# First check how 
convergence_linear %>% 
  filter(converged == FALSE) %>% 
  count(sample_size, gps_mod, outcome_interaction)


process_metrics <- function(data, model_types, figure_filename) {
  
  # collapse by mean bias at each point of exposure
  metrics_collapse <- 
    data %>%
    filter(model %in% model_types) %>%
    group_by(model, exposure, gps_mod, sample_size, outcome_interaction) %>%
    summarize(
      abs_bias = abs(mean(bias, na.rm = TRUE)),
      mse = mean(mse, na.rm = TRUE),
      rmse = sqrt(mean(mse, na.rm = TRUE)),
      .groups = 'drop'
    )
  
  labels <- if (length(model_types) == 7) {
    c("Linear", "Linear\nentropy (m = 2)", "GAM", "GAM\nentropy (m = 2)", 
      "Change", "Change\nentropy (m = 2)", "CausalGPS")
  } else {
    c("Linear\nentropy", "Linear\nentropy (m = 2)", "GAM\nentropy", "GAM\nentropy (m = 2)", 
      "Change\nentropy", "Change\nentropy (m = 2)")
  }
  
  metrics_collapse <- 
    metrics_collapse %>%
    group_by(model, gps_mod, sample_size, outcome_interaction) %>%
    summarize(
      abs_bias = mean(abs_bias, na.rm = TRUE),
      mse = mean(mse, na.rm = TRUE),
      rmse = mean(rmse, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Model = factor(model, levels = model_types, labels = labels), 
      gps_mod = factor(case_when(
        gps_mod == 1 ~ "linear", 
        gps_mod == 2 ~ "heavy tail", 
        gps_mod == 3 ~ "nonlinear", 
        gps_mod == 4 ~ "interaction"
      ), levels = c("linear", "heavy tail", "nonlinear", "interaction")), 
      out_relationship = if_else(as.logical(outcome_interaction), "interaction", "linear"), 
      out_relationship = factor(out_relationship, levels = c("linear", "interaction"))
    )
  
  # Add correct labels  
  levels(metrics_collapse$out_relationship) = c(linear = TeX("$\\mu_{linear}$"), interaction = TeX("$\\mu_{linear, int}$"))
  
  combine_plot <- bias_RMSE_combine_plot(metrics_collapse)
  ggsave(
    filename = figure_filename,
    plot = combine_plot,
    width = 10,
    height = 9,
    dpi = 400
  )
}

# Run for linear case ------------
# process_metrics(
#   data = metrics_linear, 
#   model_types = c("linear_model", "ent_linear2", "gam_model", "ent_gam2", "change_model", "ent_change2", "causal_gps_tuned"), 
#   figure_filename = paste0(figure_dir, "supp_xxx_linear_ent2.png")
# )

process_metrics(
  data = metrics_linear, 
  model_types = c("ent_linear", "ent_linear2", "ent_gam", "ent_gam2", "ent_change", "ent_change2"), 
  figure_filename = paste0(figure_dir, "supp_fig13_linear_ent_weight.pdf")
)

# Now run for sublinear case --------------------------
# process_metrics(
#   data = metrics_sublinear, 
#   model_types = c("linear_model", "ent_linear2", "gam_model", "ent_gam2", "change_model", "ent_change2", "causal_gps_tuned"), 
#   figure_filename = paste0(figure_dir, "supp_fig4_sublinear_ent2.png")
# )

process_metrics(
  data = metrics_sublinear, 
  model_types = c("ent_linear", "ent_linear2", "ent_gam", "ent_gam2", "ent_change", "ent_change2"), 
  figure_filename = paste0(figure_dir, "supp_fig14_sublinear_ent_weight.pdf")
)

# Now run threshold case
# process_metrics(
#   data = metrics_threshold, 
#   model_types = c("linear_model", "ent_linear2", "gam_model", "ent_gam2", "change_model", "ent_change2", "causal_gps_tuned"), 
#   figure_filename = paste0(figure_dir, "supp_fig5_threshold_ent2.png")
# )

process_metrics(
  data = metrics_threshold, 
  model_types = c("ent_linear", "ent_linear2", "ent_gam", "ent_gam2", "ent_change", "ent_change2"), 
  figure_filename = paste0(figure_dir, "supp_fig15_threshold_ent_weight.pdf")
)



## Supp figure 12 is produced from running 4_sensitivity_plots.R

## Figures 4-6 are produced in data_application/2_post_processing.R


# Around 50% of model 3 fails to converge for entropy balancing on second 
# moment with sample size of 200








