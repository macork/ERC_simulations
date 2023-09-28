# Script to look into data application setting
library(tidyverse)
library(latex2exp)
library(patchwork)

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
      out_relationship = if_else(as.logical(outcome_interaction), "interaction", "linear"), 
      out_relationship = factor(out_relationship, levels = c("linear", "interaction"))
    )
  
  # Add correct labels  
  levels(metrics_collapse$out_relationship) = c(linear = TeX("$\\mu_{linear}$"), interaction = TeX("$\\mu_{linear, int}$"))
  
  return(metrics_collapse)
}

# Change name of models for simplicity 
model_types = c("linear_model", "ent_linear",
                "gam_model", "ent_gam", 
                "change_model", "ent_change",
                "causal_gps_tuned")

# Load data for each model selection ------------------------------------------------
# First load the linear model example
exp_relationship = "linear"
adjust_confounder = T
time_stamp = "data_sensitivity_analysis"

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
time_stamp = "data_sensitivity_analysis"

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
time_stamp = "data_sensitivity_analysis"

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


# Now work on plots 

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


linear_plot <- bias_RMSE_combine_plot(metrics_collapse_linear)
sublinear_plot <- bias_RMSE_combine_plot(metrics_collapse_sublinear)
threshold_plot <- bias_RMSE_combine_plot(metrics_collapse_threshold)

# Combine all into one plot
metrics_collapse_linear <- metrics_collapse_linear %>% mutate(erc_relationship = "linear")
metrics_collapse_sublinear <- metrics_collapse_sublinear %>% mutate(erc_relationship = "sublinear")
metrics_collapse_threshold <- metrics_collapse_threshold %>% mutate(erc_relationship = "threshold")

metrics_collapse_all <- rbind(metrics_collapse_linear, metrics_collapse_sublinear, metrics_collapse_threshold)

metrics_collapse_all <- 
  metrics_collapse_all %>% 
  mutate(erc_relationship = factor(erc_relationship, levels = c("linear", "sublinear", "threshold")),
         out_relationship = if_else(as.logical(outcome_interaction), "interaction", "linear"), 
         out_relationship = factor(out_relationship, levels = c("linear", "interaction")))

# Add correct labels  
levels(metrics_collapse_all$out_relationship) = c(linear = TeX("$\\mu_{no-interaction}$"), interaction = TeX("$\\mu_{interaction}$"))


abs_bias_plot <- 
  metrics_collapse_all %>% 
  filter(sample_size %in% c(10000)) %>%
  ggplot() + 
  geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
  geom_point(aes(y = abs_bias, x = factor(erc_relationship), 
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
  scale_shape_manual(values=c(21, 10, 22, 12, 23, 9, 8)) +
  facet_grid( ~ out_relationship, scales = "free", labeller = label_parsed) + 
  scale_x_discrete(labels = c(linear = latex2exp::TeX("$\\mu_{linear}$"), 
                              sublinear = latex2exp::TeX("$\\mu_{sublinear}$"), 
                              threshold = latex2exp::TeX("$\\mu_{threshold}$"))) +
  labs(y = "Absolute bias", x = "")

  
  # Now plot the RMSE 
  rmse_plot <- 
    metrics_collapse_all %>% 
    filter(sample_size %in% c(10000)) %>%
    # filter(model != "causal_gps_default") %>% 
    # filter(confounder_setting == "simple") %>% 
    #filter(sample_size == 1000) %>% 
    ggplot() + 
    geom_hline(yintercept = 0, size = 0.15, linetype = "dashed") +
    geom_point(aes(y = rmse, x = factor(erc_relationship), 
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
    scale_x_discrete(labels = c(linear = latex2exp::TeX("$\\mu_{linear}$"), 
                                sublinear = latex2exp::TeX("$\\mu_{sublinear}$"), 
                                threshold = latex2exp::TeX("$\\mu_{threshold}$"))) +
    labs(y = "RMSE", x = "Mean function for outcome model")
  
  # wrap the two plots with brackets
  plot_a <- abs_bias_plot
  plot_b <- rmse_plot
  
  # arrange the two plots side by side
  combined_plots <- plot_a + plot_b & theme(legend.position = "bottom")
  combined_plots <- 
    combined_plots + 
    plot_layout(guides = "collect", ncol = 1) + 
    plot_annotation(tag_levels = 'a')
  
all_together_plot <- combined_plots

# Save plots
ggsave(
  filename = paste0(figure_dir, "supp_fig12_data_analysis_sensitivity.pdf"),
  plot = all_together_plot,
  width = 10,
  height = 9,
  dpi = 400
)


linear_plot <- bias_RMSE_combine_plot(metrics_collapse_linear)
sublinear_plot <- bias_RMSE_combine_plot(metrics_collapse_sublinear)
threshold_plot <- bias_RMSE_combine_plot(metrics_collapse_threshold)

