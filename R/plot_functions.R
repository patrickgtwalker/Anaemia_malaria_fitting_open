# Function to process MCMC output and generate data frames for plotting
summarise_model <- function(mcmc, site_inf_history, n_samples=100,gestage_min=80,gestage_max=200) {
  
  # Determine the total number of samples available
  max_list <- (mcmc$parameters$samples * mcmc$parameters$chains)
  # Randomly sample n_samples from the total number of samples
  n_list <- sample(1:max_list, n_samples, replace = FALSE) 
  
  # Extract the MCMC output
  output <- mcmc$output
  
  # Generate malaria splines for each sample
  mal_splines_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("mal_knot_%s", 1:3)) %>%
    apply(MARGIN = 1, FUN = function(y) cubic_spline(x = seq(gestage_min, gestage_max, l = 3), y = y, 
                                                     x_pred = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)))) %>%
    as.data.frame() %>%
    mutate(x = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)))
  
  # Generate gestational age (GA) splines for each sample
  GA_splines_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("y_knot_%s", 1:3)) %>%
    apply(MARGIN = 1, FUN = function(y) cubic_spline(x = seq(gestage_min, gestage_max, l = 3), y = y, 
                                                     x_pred = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)))) %>%
    as.data.frame() %>%
    mutate(x = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)))
  
  # Create an array to store G_non_infect values for each gravidity category (G1 to G6)
  G_non_infect_array <- array(0, dim = c(max_list, 6))  # 6 columns for G1 to G6
  
  # Fill the array with G_non_infect values for each gravidity category
  for (i in 2:6) {
    G_non_infect_array[, i] <- output %>%
      filter(phase == "sampling") %>%
      select(sprintf("G%s_non_infect", i)) %>%
      unlist()  # Convert to a vector and store in the array
  }
  
  # Extract the log odds of malaria prevalence in G1
  log_odds_malaria_prevalence_G1_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("log_odds_malaria_prevalence_G1")) %>%
    as.data.frame()
  
  # Extract the second malaria knot value
  mal_knot2 <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("mal_knot_2")) %>%
    as.data.frame()
  
  # Extract the Hill function scale parameter
  scale_hill_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("scale_hill")) %>%
    as.data.frame()
  
  # Extract the Hill function shape parameter
  shape_hill_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("shape_hill")) %>%
    as.data.frame()
  
  # Initialize empty lists to store data frames for each draw, gravidity, and malaria status
  HB_mean_list <- list()
  mal_splines_list <- list()
  mal_immunity_list <- list()
  
  # Loop over the draws in n_list
  for (j in n_list) {
    # Generate the immunity vector using the Hill function
    imm_vect <- hill_func(0:20, shape = shape_hill_out$shape[j], scale = scale_hill_out$scale[j])
    
    # Calculate gravidity impact based on the first gravidity's malaria prevalence
    grav_impact_df <- get_weight_impact_df(inf_history = site_inf_history,
                                           primi_prev = plogis(log_odds_malaria_prevalence_G1_out$log_odds_malaria_prevalence_G1[j]),
                                           imm_vect = imm_vect)
    
    # Create a data frame for malaria splines for this draw
    mal_splines_df <- data.frame(
      mal_effect = -mal_splines_out[, j],
      gestage = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)),
      draw = j
    )
    
    # Loop over gravidities (1 to 6)
    for (g in 1:6) {
      # Calculate malaria effect by gravidity
      mal_immunity <- data.frame(
        mal_effect = -mal_knot2$mal_knot_2[j] * grav_impact_df[g],
        draw = j,  # Store the draw index
        gravidity = g
      )
      
      # Calculate HB_mean with malaria for each gravidity
      HB_mean_mal <- data.frame(
        HB_mean = GA_splines_out[, j] + mal_splines_out[, j] * grav_impact_df[g] + G_non_infect_array[j, g],
        gestage = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)),
        draw = j,  # Store the draw index
        gravidity = g,
        malaria = 1  # Indicate malaria status
      )
      
      # Calculate HB_mean without malaria for each gravidity
      HB_mean_no_mal <- data.frame(
        HB_mean = GA_splines_out[, j] + G_non_infect_array[j, g],
        gestage = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)),
        draw = j,  # Store the draw index
        gravidity = g,
        malaria = 0  # Indicate no malaria status
      )
      
      # Append the data frames to the lists
      mal_immunity_list[[length(mal_immunity_list) + 1]] <- mal_immunity
      mal_splines_list[[length(mal_splines_list) + 1]] <- mal_splines_df
      HB_mean_list[[length(HB_mean_list) + 1]] <- HB_mean_mal
      HB_mean_list[[length(HB_mean_list) + 1]] <- HB_mean_no_mal
    }
  }
  
  # Combine all the data frames into one using bind_rows
  HB_mean_df <- bind_rows(HB_mean_list)
  mal_splines_df <- bind_rows(mal_splines_list)
  mal_immunity_df <- bind_rows(mal_immunity_list)
  
  # Rename the value columns
  colnames(HB_mean_df)[1] <- "HB_mean"
  colnames(mal_splines_df)[1] <- "mal_effect"
  colnames(mal_immunity_df)[1] <- "mal_effect"
  
  # Return the combined data frames as a list
  return(list(
    mal_splines_df = mal_splines_df,
    HB_mean_df = HB_mean_df,
    mal_immunity_df = mal_immunity_df
  ))
}


generate_combined_plot <- function(model, fitted_data,true_data,site_inf_history, n_samples=1000,gestage_min=80,gestage_max=200){

# Fit the non-censored data
fit_summary <- summarise_model(mcmc=model,site_inf_history=site_inf_history, n_samples=n_samples)
# First plot: HB Mean vs. Gestational Age by Gravidity and Malaria Status
plot1_noncensored <- ggplot(fit_summary$HB_mean_df) +
  geom_point(data = fitted_data, aes(x = gestage, y = hb_level, color = factor(malaria)), alpha = 0.05) +
  geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria)), alpha = 0.1) +
  stat_smooth(data = true_data$data_df, aes(x = gestage, y = hb_level, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dotted", linewidth = 2) +
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +
  facet_wrap(~ gravidity) +
  theme_minimal()

# Second plot: Malaria attributable reduction in HB over gestational age
plot2_noncensored <- ggplot(fit_summary$mal_splines_df) +
  geom_line(aes(x = gestage, y = mal_effect, group = draw), alpha = 0.05,color = "#6B5072") +
  geom_line(data = true_data$simmed_mal_spline, aes(x = gestage, y = mal_effect), linewidth = 3,color = "#D08C60") +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()

# Third plot: Malaria attributable reduction in HB over gravidity
plot3_noncensored <- ggplot(fit_summary$mal_immunity_df) +
  geom_line(aes    (x = gravidity, y = mal_effect, group = draw), alpha = 0.05,color = "#6B5072") +
  geom_line(data = true_data$simmed_mal_immunity, aes(x = gravidity, y = mal_effect), linewidth = 3,color = "#D08C60") +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()

# Combine the three plots into a dashboard
return(wrap_elements(plot1_noncensored) + 
  (wrap_elements(plot2_noncensored) / wrap_elements(plot3_noncensored)) + 
  plot_layout(widths = c(2, 1))
)
}
