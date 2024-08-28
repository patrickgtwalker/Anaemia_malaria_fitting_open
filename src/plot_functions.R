library("patchwork")
mcmc_out_new <- function(mcmc,n_samples) {
  
  max_list <- (mcmc$parameters$samples * mcmc$parameters$chains)
  n_list <- sample(1:max_list, n_samples, replace=F) #if only sample for n_list once, it stays the same.
  output<-mcmc$output
  mal_splines_out<- output %>%
    filter(phase=="sampling") %>%
    select(sprintf("mal_knot_%s", 1:3)) %>%
    apply(MARGIN = 1, FUN = function(y) cubic_spline(x = seq(gestage_min,gestage_max,l=3), y=y, 
                                                     x_pred = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1))))%>%
    as.data.frame() %>%
    mutate(x = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1)))
  
  GA_splines_out<- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("y_knot_%s", 1:3)) %>%
    apply(MARGIN = 1, FUN = function(y) cubic_spline(x = seq(gestage_min,gestage_max,l=3), y=y, 
                                                     x_pred = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1))))%>%
    as.data.frame() %>%
    mutate(x = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1)))
  
  G_non_infect_array <- array(0, dim = c(max_list, 6))  # 6 columns for G1 to G6
  
  # Fill the array with G_non_infect values
  for (i in 2:6) {
    G_non_infect_array[,i] <- output %>%
      filter(phase == "sampling") %>%
      select(sprintf("G%s_non_infect", i)) %>%
      unlist()  # Convert to a vector and store in the array
  }
  log_odds_malaria_prevalence_G1_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("log_odds_malaria_prevalence_G1")) %>%
    as.data.frame()
  
  mal_knot2 <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("mal_knot_2")) %>%
    as.data.frame()
 
   scale_hill_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("scale_hill")) %>%
    as.data.frame()
  
  shape_hill_out <- output %>%
    filter(phase == "sampling") %>%
    select(sprintf("shape_hill")) %>%
    as.data.frame()
  
  
  
  # Initialize an empty list to store data frames for each draw, gravidity, and malaria status
  HB_mean_list <- list()
  mal_splines_list <- list()
  mal_immunity_list <- list()
  # Loop over the draws in n_list
  for(j in n_list){
    # Generate the immunity vector
    imm_vect <- hill_func(0:20, shape = shape_hill_out$shape[j], scale = scale_hill_out$scale[j])
    
    # Calculate gravidity impact based on the first gravidity's malaria prevalence
    grav_impact_df <- get_weight_impact_df(inf_history = burkina_inf_history,
                                           primi_prev = plogis(log_odds_malaria_prevalence_G1_out$log_odds_malaria_prevalence_G1[j]),
                                           imm_vect = imm_vect)
    mal_splines_df<-data.frame(mal_effect =-mal_splines_out[,j],
               gestage = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)),
               draw = j)
    # Loop over gravidities (assuming 1 to 6 for gravidity categories)
    for(g in 1:6){
      # Calculate HB_mean with malaria for each gravidity
      mal_immunity<-data.frame(mal_effect =-mal_knot2$mal_knot_2[j] * grav_impact_df[g],
                               draw = j,  # Store the draw index
                               gravidity = g
                               )
      HB_mean_mal <- data.frame(HB_mean = GA_splines_out[,j] + mal_splines_out[,j] * grav_impact_df[g] + G_non_infect_array[j, g],
                                gestage = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)),
                                draw = j,  # Store the draw index
                                gravidity = g,
                                malaria = 1)  # Indicate malaria status
      
      # Calculate HB_mean without malaria for each gravidity
      HB_mean_no_mal <- data.frame(HB_mean = GA_splines_out[,j] + G_non_infect_array[j, g],
                                   gestage = seq(gestage_min, gestage_max, l = ((gestage_max - gestage_min) + 1)),
                                   draw = j,  # Store the draw index
                                   gravidity = g,
                                   malaria = 0)  # Indicate no malaria status
      
      # Append the data frames to the list
      mal_immunity_list[[length(mal_immunity_list) + 1]] <- mal_immunity
      mal_splines_list[[length(mal_splines_list) + 1]] <- mal_splines_df
      HB_mean_list[[length(HB_mean_list) + 1]] <- HB_mean_mal
      HB_mean_list[[length(HB_mean_list) + 1]] <- HB_mean_no_mal
    }
  }
  
  # Combine all the data frames into one using bind_rows
  HB_mean_df <- bind_rows(HB_mean_list)
  mal_splines_df<-bind_rows(mal_splines_list)
  mal_immunity_df<-bind_rows(mal_immunity_list)
  # Rename the value column (assuming the column for HB_mean_mal is V1)
  colnames(HB_mean_df)[1] <- "HB_mean"
  colnames(mal_splines_df)[1]<-"mal_effect"
  colnames(mal_immunity_df)[1]<-"mal_effect"
  
  return(list(mal_splines_df=mal_splines_df,
            HB_mean_df=HB_mean_df,
            mal_immunity_df=mal_immunity_df))
}

noncensored_fit<-mcmc_out_new(run_noncensored_data,100)
plot1_noncensored <- ggplot(noncensored_fit$HB_mean_df ) +
  geom_point(data = synthetic_data_df, aes(x = gestage, y = hb_level, color = factor(malaria)), alpha = 0.05) +
  geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria)),alpha=0.01) +
  stat_smooth(data = synthetic_data_df, aes(x = gestage, y = hb_level, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dotted", size = 2) +
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +
  facet_wrap(~ gravidity) +
  theme_minimal()


# Second plot: Malaria attributable reduction in HB over gestational age
plot2_noncensored <- ggplot(noncensored_fit$mal_splines_df) +
  geom_line(aes(x = gestage, y = mal_effect, group = draw), alpha = 0.1) +
  geom_line(data = simmed_mal_spline, aes(x = gestage, y = mal_effect), linewidth = 3) +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()

# Third plot: Malaria attributable reduction in HB over gravidity
plot3_noncensored <- ggplot(noncensored_fit$mal_immunity_df) +
  geom_line(aes(x = gravidity, y = mal_effect, group = draw), alpha = 0.1) +
  geom_line(data = simmed_mal_immunity, aes(x = gravidity, y = mal_effect), linewidth = 3) +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()

combined_noncensored_plot <- wrap_elements(plot1_noncensored) + 
  (wrap_elements(plot2_noncensored) / wrap_elements(plot3_noncensored)) + 
  plot_layout(widths = c(2, 1))

combined_noncensored_plot

incorrect_censored_fit<-mcmc_out_new(run_censored_data_incorrect,100)
plot1_incorrect_censored <- ggplot(incorrect_censored_fit$HB_mean_df ) +
  geom_point(data = synthetic_data_df_censor, aes(x = gestage, y = hb_level, color = factor(malaria)), alpha = 0.05) +
  geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria)),alpha=0.01) +
  stat_smooth(data = synthetic_data_df, aes(x = gestage, y = hb_level, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dotted", size = 2) +
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +
  facet_wrap(~ gravidity) +
  theme_minimal()


# Second plot: Malaria attributable reduction in HB over gestational age
plot2_incorrect_censored <- ggplot(incorrect_censored_fit$mal_splines_df) +
  geom_line(aes(x = gestage, y = mal_effect, group = draw), alpha = 0.1) +
  geom_line(data = simmed_mal_spline, aes(x = gestage, y = mal_effect), linewidth = 3) +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()

# Third plot: Malaria attributable reduction in HB over gravidity
plot3_incorrect_censored <- ggplot(incorrect_censored_fit$mal_immunity_df) +
  geom_line(aes(x = gravidity, y = mal_effect, group = draw), alpha = 0.1) +
  geom_line(data = simmed_mal_immunity, aes(x = gravidity, y = mal_effect), linewidth = 3) +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()

combined_incorrect_censored_plot <- wrap_elements(plot1_incorrect_censored ) + 
  (wrap_elements(plot2_incorrect_censored ) / wrap_elements(plot3_incorrect_censored )) + 
  plot_layout(widths = c(2, 1))

combined_incorrect_censored_plot


correct_censored_fit<-mcmc_out_new(run_censored_data_correct,100)
plot1_correct_censored <- ggplot(correct_censored_fit$HB_mean_df ) +
  geom_point(data = synthetic_data_df_censor, aes(x = gestage, y = hb_level, color = factor(malaria)), alpha = 0.05) +
  geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria)),alpha=0.01) +
  stat_smooth(data = synthetic_data_df, aes(x = gestage, y = hb_level, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dotted", size = 2) +
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +
  facet_wrap(~ gravidity) +
  theme_minimal()


# Second plot: Malaria attributable reduction in HB over gestational age
plot2_correct_censored <- ggplot(correct_censored_fit$mal_splines_df) +
  geom_line(aes(x = gestage, y = mal_effect, group = draw), alpha = 0.1) +
  geom_line(data = simmed_mal_spline, aes(x = gestage, y = mal_effect), linewidth = 3) +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()

# Third plot: Malaria attributable reduction in HB over gravidity
plot3_correct_censored <- ggplot(correct_censored_fit$mal_immunity_df) +
  geom_line(aes(x = gravidity, y = mal_effect, group = draw), alpha = 0.1) +
  geom_line(data = simmed_mal_immunity, aes(x = gravidity, y = mal_effect), linewidth = 3) +
  ylab("Malaria attributable reduction in HB (g/dL)") +
  theme_minimal()


combined_correct_censored_plot <- wrap_elements(plot1_correct_censored ) + 
  (wrap_elements(plot2_correct_censored ) / wrap_elements(plot3_correct_censored )) + 
  plot_layout(widths = c(2, 1))

combined_correct_censored_plot

