# Function to fit CausalGPS by each year of the dataset
# nrounds, max depth, eta are parameters for xgboost and fittings GPS
# delta specifies bin length for expsoure, another hyperparameters
# Return_Corr_only only returns the correlations if looking to train hyperparameters only
# data ent is the data used for creating gps 
# Confounders is list of confounder names 
fit_causalGPS_by_year <- function(nrounds,
                                  max_depth,
                                  eta,
                                  delta,
                                  return_corr_only = T,
                                  nthread = 1,
                                  trim = T,
                                  data_ent, 
                                  confounders) {
  
  
  # Add ID column for CausalGPS package
  data_ent <- mutate(data_ent, id = row_number())
  data_confounders <- data.frame(data_ent %>% select(all_of(c(confounders, "id"))) %>% select(-year))
  
  # Estimate the GPS
  zip_year_gps <- 
    estimate_gps(w = data.frame(dplyr::select(data_ent, pm25, id)),
                 c = data_confounders,
                 ci_appr = "matching",
                 pred_model = "sl",
                 gps_density = "normal",
                 use_cov_transform = FALSE,
                 transformers = list("pow2", "pow3", "abs", "scale"),
                 trim_quantiles = c(0, 1),
                 optimized_compile = TRUE,
                 sl_lib = c("m_xgboost"),
                 params = list(
                   "xgb_nrounds" = nrounds,
                   "xgb_max_depth" = max_depth,
                   "xgb_eta" = eta,
                   "xgb_min_child_weight" = 1
                 ),
                 covar_bl_method = "absolute",
                 covar_bl_trs = 0.1,
                 covar_bl_trs_type = "mean",
                 max_attempt = 1,
                 matching_fun = "matching_l1",
                 delta_n = delta,
                 scale = 1,
                 nthread = nthread)
  
  # Extract data to then split by year with GPS attached
  # Group by year
  zip_year_gps_plus_params <- 
    zip_year_gps$dataset %>%
    left_join(data_ent)
  
  # Split the data frame into a list of data frames, one for each year
  zip_list_by_year <- 
    zip_year_gps_plus_params %>%
    split(.$year)
  
  # now match within each year
  matched_data <- 
    rbindlist(lapply(zip_list_by_year, function(z){
      
      # Create dataset so CausalGPS can match
      dataset_gps <- list()
      class(dataset_gps) <- "cgps_gps"
      
      dataset_gps$dataset <- 
        z %>%
        select(id, pm25, gps, e_gps_pred, e_gps_std_pred, w_resid, zip, year, all_of(confounders)) %>% 
        as.data.frame()
      dataset_gps$gps_mx <- zip_year_gps$gps_mx
      dataset_gps$w_mx <- zip_year_gps$w_mx
      
      # now match dataset
      matched_pop <- compile_pseudo_pop(data_obj = dataset_gps,
                                        ci_appr = "matching",
                                        gps_density = "normal",
                                        exposure_col_name = c("pm25"),
                                        bin_seq = NULL,
                                        nthread = nthread,
                                        dist_measure = "l1",
                                        covar_bl_method = "absolute",
                                        covar_bl_trs = 0.1,
                                        covar_bl_trs_type = "maximal",
                                        delta_n = delta,
                                        scale = 1)
      
      return(matched_pop)
    }))
  
  # Truncate matching weights to 99th percentile (if specified)
  if (trim) {
    matched_data[counter_weight > quantile(matched_data$counter_weight, 0.99),
                 counter_weight := quantile(matched_data$counter_weight, 0.99)] 
  }
  
  # Now generate covariate balance tab
  balance_table <- 
    bal.tab(reformulate(confounders, response = "pm25"),
            data = matched_data,
            weights = matched_data[["counter_weight"]],
            method = "weighting",
            stats = c("cor"),
            un = T,
            continuous = "std",
            s.d.denom = "weighted",
            abs = T,
            thresholds = c(cor = .1), poly = 1)
  
  mean_post_cor <- mean(balance_table$Balance$Corr.Adj)
  
  # Save absolute correlation
  corr_results <- data.table(nrounds = nrounds, depth = max_depth,
                             eta = eta, delta = delta,
                             corr = mean_post_cor)
  
  # If just fitting correlation return here
  if (return_corr_only) {
    return(corr_results)
  }
  # Return psuedo pop, balance table, correlation resluts
  results <- list(pseudo_pop = matched_data,
                  balance_table = balance_table,
                  corr_results = corr_results)
  return(results)
}
