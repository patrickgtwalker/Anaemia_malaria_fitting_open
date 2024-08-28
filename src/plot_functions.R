
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
  

uncensored_fit<-mcmc_out_new(run_noncensored_data,100)
ggplot(uncensored_fit$HB_mean_df ) +
  geom_point(data=synthetic_data_df,aes(x=gestage,y=hb_level, color = factor(malaria)),alpha=0.05)+
  geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria))) +
  stat_smooth(data=synthetic_data_df,aes(x=gestage,y=hb_level, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dashed", size = 1)+
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +
  facet_wrap(~ gravidity) +  # Facet wrap by gravidity
  theme_minimal()


incorrect_censored_fit<-mcmc_out_new(run_censored_data_incorrect,100)
correct_censored_fit<-mcmc_out_new(run_censored_data_correct,100)

ggplot(correct_fit$mal_splines_df) +
  geom_line(aes(x = gestage, y = mal_effect, group = draw),alpha=0.1)+
  geom_line(data=simmed_mal_spline,aes(x = gestage, y = mal_effect),linewidth=3)
  
ggplot(correct_fit$mal_immunity_df) +
  geom_line(aes(x = gravidity, y = mal_effect, group = draw),alpha=0.1)+
  geom_line(data=simmed_mal_immunity,aes(x = gravidity, y = mal_effect),linewidth=3)

ggplot(incorrect_fit$mal_immunity_df) +
  geom_line(aes(x = gravidity, y = mal_effect, group = draw),alpha=0.1)+
  geom_line(data=simmed_mal_immunity,aes(x = gravidity, y = mal_effect),linewidth=3)+ylab("malari")

ggplot(incorrect_fit$mal_splines_df) +
  geom_line(aes(x = gestage, y = mal_effect, group = draw),alpha=0.1)+
  geom_line(data=simmed_mal_spline,aes(x = gestage, y = mal_effect),linewidth=3)


ggplot(incorrect_fit$HB_mean_df ) +
  geom_point(data=synthetic_data_df,aes(x=gestage,y=hb_level, color = factor(malaria)),alpha=0.05)+
    geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria))) +
  stat_smooth(data=synthetic_data_df,aes(x=gestage,y=hb_level, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dashed", size = 1)+
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +

  facet_wrap(~ gravidity) +  # Facet wrap by gravidity
  theme_minimal()



ggplot(correct_fit$HB_mean_df ) +
  geom_point(data=synthetic_data_df,aes(x=gestage,y=hb_level, color = factor(malaria)),alpha=0.05)+
  geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria))) +
  stat_smooth(data=synthetic_data_df,aes(x=gestage,y=hb_level, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dashed", size = 1)+
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +
  
  facet_wrap(~ gravidity) +  # Facet wrap by gravidity
  theme_minimal()


length(check$GA_splines_out[,1])

quick_run$output%>%
  filter(phase=="sampling") %>%
  select(sprintf("mal_knot_%s", 1:3)) %>%
  apply(MARGIN = 1, FUN = function(y) cubic_spline(x = seq(gestage_min,gestage_max,l=3), y=y, 
                                                   x_pred = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1))))%>%
  as.data.frame() %>%
  mutate(x = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1)))


mcmc_out_new(quick_run$output)
 
  GA_splines_out<- mcmc %>%
    filter(phase == "sampling") %>%
    select(sprintf("y_knot_%s", 1:3)) %>%
    apply(MARGIN = 1, FUN = function(y) cubic_spline(x = seq(min,max,l=3), y=y, 
                                                     x_pred = seq(min,max,l=(2*(max-min)+1)))) %>%
    as.data.frame() %>%
    mutate(x =  seq(min,max,l=(2*(max-min)+1)))

 
  HB_sigma_out <- mcmc %>%
    filter(phase == "sampling") %>%
    select(sprintf("HB_sigma")) %>%
    as.data.frame()
  
  
  for (i in 2:6) {
    assign(paste0("G", i, "_non_infect_out"), 
           mcmc %>%
             filter(phase == "sampling") %>%
             select(sprintf("G%s_non_infect", i)) %>%
             as.data.frame())
  }
  
  scale_hill_out <- mcmc %>%
    filter(phase == "sampling") %>%
    select(sprintf("scale_hill")) %>%
    as.data.frame()
  
  shape_hill_out <- mcmc %>%
    filter(phase == "sampling") %>%
    select(sprintf("shape_hill")) %>%
    as.data.frame()
  
  
  return(list(GA_splines_out = GA_splines_out, mal_splines_out = mal_splines_out,
              site_out_1 = site_out_1, site_out_2 = site_out_2, site_out_3 = site_out_3,
              site_out_4 = site_out_4, site_out_5 = site_out_5, site_out_6 = site_out_6,
              site_out_7 = site_out_7, site_out_8 = site_out_8, 
              log_odds_PG_prev_1_out = log_odds_PG_prev_1_out, log_odds_PG_prev_2_out = log_odds_PG_prev_2_out,
              log_odds_PG_prev_3_out = log_odds_PG_prev_3_out, log_odds_PG_prev_4_out = log_odds_PG_prev_4_out,
              log_odds_PG_prev_5_out = log_odds_PG_prev_5_out, log_odds_PG_prev_6_out = log_odds_PG_prev_6_out,
              log_odds_PG_prev_7_out = log_odds_PG_prev_7_out, log_odds_malaria_prevalence_G1_wkenya_out = log_odds_malaria_prevalence_G1_wkenya_out,
              log_odds_malaria_prevalence_G1_malawi_out = log_odds_malaria_prevalence_G1_malawi_out, HB_sigma_out = HB_sigma_out, 
              G2_non_infect_out = G2_non_infect_out, G3_non_infect_out = G3_non_infect_out,
              G4_non_infect_out = G4_non_infect_out, G5_non_infect_out = G5_non_infect_out,
              G6_non_infect_out = G6_non_infect_out,
              scale_hill_out = scale_hill_out, shape_hill_out = shape_hill_out))
} #