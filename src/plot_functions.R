n_samples <- 100 # how many samples do I want
max_list <- 2*(mcmc_G1_G2_G3_G4_G5_G6_1e41e5_model_1$parameters$samples * mcmc_G1_G2_G3_G4_G5_G6_1e41e5_model_1$parameters$chains)
n_list <- sample(1:max_list, n_samples, replace=F) #if only sample for n_list once, it stays the same.


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
  
  # Loop over the draws in n_list
  for(j in n_list){
    # Generate the immunity vector
    imm_vect <- hill_func(0:20, shape = shape_hill_out$shape[j], scale = scale_hill_out$scale[j])
    
    # Calculate gravidity impact based on the first gravidity's malaria prevalence
    grav_impact_df <- get_weight_impact_df(inf_history = burkina_inf_history,
                                           primi_prev = plogis(log_odds_malaria_prevalence_G1_out$log_odds_malaria_prevalence_G1[j]),
                                           imm_vect = imm_vect)
    
    # Loop over gravidities (assuming 1 to 6 for gravidity categories)
    for(g in 1:6){
      # Calculate HB_mean with malaria for each gravidity
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
      HB_mean_list[[length(HB_mean_list) + 1]] <- HB_mean_mal
      HB_mean_list[[length(HB_mean_list) + 1]] <- HB_mean_no_mal
    }
  }
  
  # Combine all the data frames into one using bind_rows
  HB_mean_df <- bind_rows(HB_mean_list)
  
  
  # Rename the value column (assuming the column for HB_mean_mal is V1)
  colnames(HB_mean_df)[1] <- "HB_mean"
  

  

  #  HB_mean_no_mal_G1 <- (rep(site_out_x[j,], rep_number) + mcmc_output$GA_splines_out[,j]) 
 # HB_mean_mal_G1 <- as.data.frame(HB_mean_mal_G1) %>% mutate(x=seq(min,max,l=(2*(max-min)+1))) 
  #HB_mean_no_mal_G1 <- as.data.frame(HB_mean_no_mal_G1) %>% mutate(x=seq(min,max,l=(2*(max-min)+1)))
  #HB_mean_mal_G1 <- cbind(HB_mean_mal_G1, rep(j, rep_number), rep(1, rep_number), rep(country_name, rep_number), rep(1, rep_number))
  #HB_mean_no_mal_G1 <- cbind(HB_mean_no_mal_G1, rep(j, rep_number), rep(1, rep_number), rep(country_name, rep_number), rep(0, rep_number))
  #malaria_output_G1 = rbind(malaria_output_G1, HB_mean_mal_G1)
  #no_malaria_output_G1 = rbind(no_malaria_output_G1, HB_mean_no_mal_G1)
  
  
  return(list(mal_splines_out=mal_splines_out,
              GA_splines_out=GA_splines_out,
              HB_mean_df=HB_mean_df))
}
  

check<-mcmc_out_new(run_censored_data_incorrect,100)

ggplot(check$HB_mean_df ) +
  geom_point(data=synthetic_data_df,aes(x=gestage,y=hb_level, color = factor(malaria)),alpha=0.1)+
    geom_line(aes(x = gestage, y = HB_mean, group = interaction(malaria, draw), color = factor(malaria))) +
  stat_smooth(aes(x = gestage, y = HB_mean, color = factor(malaria)), method = "loess", se = FALSE, linetype = "dashed", size = 1)+
  labs(title = "HB Mean vs. Gestational Age by Gravidity and Malaria Status",
       x = "Gestational Age",
       y = "HB Mean",
       color = "Malaria Status") +

  facet_wrap(~ gravidity, scales = "free_y") +  # Facet wrap by gravidity
  theme_minimal()


check$HB_mean_mal_G1_combined
ggplot(check$HB_mean_df) +
  geom_line(aes(x = gestage, y = HB_mean, group = draw)) +
  labs(title = "HB Mean Malaria Effect across Gestational Age",
       x = "Gestational Age",
       y = "HB Mean") +
  theme_minimal()+ylim(0,20)+
  stat_smooth(aes(x = gestage, y = HB_mean), method = "loess", se = FALSE, linetype = "dashed", size = 1)+
  geom_point(data=synthetic_data_df%>%filter(gravidity==1,malaria==1),aes(x=gestage,y=hb_level))

head(synthetic_data_df)

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