# Calculate the distribution of gravidity categories based on relative risks (RRs)
get_dist_grav <- function(RRs) {
  numerator <- c(1, RRs)
  denom <- sum(numerator)
  return(numerator / denom)
}
# Function to compute immunity vector using the Hill function
hill_func <- function(i, shape, scale) {
  return(1 / (1 + (i / scale)^shape))
}

generate_synthetic_data <- function(
    gestage_min = 80, 
    gestage_max = 200, 
    base_sample_size = 30000, 
    grav_cats = 6, 
    GA_shape_1 = 1, 
    GA_shape_2 = 1, 
    log_RR = c(-0.1, -0.2, -0.3, -0.4, 0.2), 
    Primigravid_prevalence = 0.5, 
    HB_sigma = 1.8, 
    y_knots = c(12, 10.5, 9.5), 
    mal_knots = c(-0.5, -1.5, -1.25), 
    hill_shape = 3, 
    hill_scale = 0.5, 
    G_non_infect_params = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), 
    site_inf_history
) {

  # Calculate the relative distribution of sample sizes for each gravidity category using log_RR
  grav_dist <- get_dist_grav(exp(log_RR))
  
  # Generate sample sizes for each gravidity category using a multinomial distribution
  sample_sizes <- as.vector(rmultinom(1, base_sample_size, grav_dist))
  
  # Calculate base log odds and adjust for each gravidity category
  log_odds_base <- log(Primigravid_prevalence / (1 - Primigravid_prevalence))
  log_odds <- log_odds_base + c(0, -0.2, -0.4, -0.6, -0.8, -1.0)
  
  # Convert log odds to prevalence
  prevalence_by_grav <- exp(log_odds) / (1 + exp(log_odds))
  
  # Generate y_spline and mal_spline values based on the cubic spline function
  y_spline <- cubic_spline(x = seq(gestage_min, gestage_max, length.out = 3), y = y_knots, x_pred = seq(gestage_min, gestage_max, length.out = (gestage_max - gestage_min) + 1))
  mal_spline <- cubic_spline(x = seq(gestage_min, gestage_max, length.out = 3), y = mal_knots, x_pred = seq(gestage_min, gestage_max, length.out = (gestage_max - gestage_min) + 1))
  
  # Generate immunity vector based on the Hill function
  imm_vect <- hill_func(0:20, hill_shape, hill_scale)
  
  # Calculate the immunity effects for each gravidity category
  imm_effects <- get_weight_impact_df(site_inf_history, Primigravid_prevalence, imm_vect)
  
  # Initialize empty lists to hold the synthetic data
  synthetic_data <- list()
  
  # Loop through each gravidity category to generate synthetic data
  for (i in 1:grav_cats) {
    # Generate gestational ages for each individual
    adj_gestage <- floor(rbeta(n = sample_sizes[i], shape1 = GA_shape_1, shape2 = GA_shape_2) * (gestage_max - gestage_min)) + 1
    gestage <- adj_gestage + gestage_min - 1
    
    # Generate malaria status (1 = malaria, 0 = no malaria) based on the prevalence
    malaria_status <- rbinom(sample_sizes[i], 1, prevalence_by_grav[i])
    
    # Access spline values based on rounded gestational age
    y_spline_values <- y_spline[adj_gestage]
    mal_spline_values <- mal_spline[adj_gestage]
    
    # Generate hemoglobin levels or other outcomes based on G_non_infect and malaria status
    hb_levels <- rnorm(sample_sizes[i], mean = G_non_infect_params[i] + y_spline_values + mal_spline_values * malaria_status * imm_effects[i], sd = HB_sigma)
    
    # Create a data frame for this gravidity category
    synthetic_data[[i]] <- data.frame(
      gravidity = i,
      gestage = gestage,
      adj_gestage = adj_gestage,
      malaria = malaria_status,
      hb_level = hb_levels
    )
  }
  simmed_mal_spline <- data.frame(
    gestage = seq(gestage_min, gestage_max, length.out = (gestage_max - gestage_min) + 1),
    mal_effect = -mal_spline  # Negate the malaria effect
  )
  # Combine all gravidity categories into one data frame
  synthetic_data_df <- do.call(rbind, synthetic_data)
  simmed_mal_immunity <- data.frame(
    gravidity = 1:6,
    mal_effect = -mal_knots[2] * imm_effects  # Negate and scale malaria effects by immunity
  )   
  return(
    list(
      data_df=synthetic_data_df,
      simmed_mal_spline=simmed_mal_spline,
      simmed_mal_immunity=simmed_mal_immunity
         )
  )
}
generate_country_data_list_noncensor <- function(data_df) {
  # Calculate data on overall risk in primigravidae for computational efficiency
  PG_prev <- data_df %>%
    filter(gravidity == 1) %>%
    summarise(positive = sum(malaria), total = n())
  
  # Create a list of the different datasets for non-censored data
  country_data_list_noncensor <- list(
    country_df = data_df,
    prev_data = PG_prev
  )
  
  return(country_data_list_noncensor)
}

generate_country_data_list_censor <- function(data_df, cutoff = 7) {
  # Mimic censored data where women with HB < 7 were excluded
  data_df_censor <- data_df %>% filter(hb_level > cutoff)
  
  # Provide data on the number of women censored
  censored_number <- as.numeric(data_df %>% filter(hb_level <= cutoff) %>% summarise(n()))
  uncensored_number <- nrow(data_df_censor)
  censored_data <- list(
    censored_number = censored_number,
    total_number = censored_number + uncensored_number
  )
  
  # Pre-calculate the observed gestational age distribution in uncensored women
  uncensored_dist <- data_df_censor %>%
    group_by(adj_gestage, malaria, gravidity, .drop = FALSE) %>%
    summarise(count = n())
  
  # Index based on gravidity and malaria status to align with function in prob_censor.R
  uncensored_dist$index <- with(uncensored_dist, (adj_gestage - 1) * (max(uncensored_dist$malaria) + 1) * max(uncensored_dist$gravidity) + malaria * max(uncensored_dist$gravidity) + gravidity)
  
  # Observed distribution by gravidity and malaria status
  grav_mal_counts <- data_df_censor %>% 
    group_by(gravidity, malaria) %>%
    summarise(total = n(), prop = n() / uncensored_number)
  
  PG_prev <- data_df %>%
    filter(gravidity == 1) %>%
    summarise(positive = sum(malaria), total = n())
  
  # Create a list of the different datasets for censored data
  country_data_list_censor <- list(
    country_df = data_df_censor,
    prev_data = PG_prev,
    grav_mal_counts = grav_mal_counts,
    counts_df = uncensored_dist,
    censored_data = censored_data
  )
  
  return(country_data_list_censor)
}
