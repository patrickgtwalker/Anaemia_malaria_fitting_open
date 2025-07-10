# Function to process the censoring block and calculate the likelihood
process_censor_block <- function(counts_df, censored_data, min_day, max_day, cutoff, max_malaria, site_FE, y_spline, mal_spline, HB_sigma, G_non_infect, grav_impact_df, grav_cats, shape_1_GA, shape_2_GA, malaria_prevalence, RRs, grav_mal_numbers) {
  
  # Compute the probability of being censored based on the input parameters
  prob_cens_with_params <- prob_cens(
    counts_df = counts_df, min_day = min_day, max_day = max_day, cutoff = cutoff, max_malaria = max_malaria, 
    site_FE = site_FE, y_spline = y_spline, mal_spline = mal_spline, HB_sigma = HB_sigma, G_non_infect = G_non_infect, grav_impact_df = grav_impact_df,
    grav_cats = grav_cats, shape_1_GA = shape_1_GA, shape_2_GA = shape_2_GA
  )
  
  # Retrieve the likelihood of observing gestational age (GA) given malaria status and gravidity
  GA_likelihood <- prob_cens_with_params$GA_likelihood
  
  # Calculate gravidity probabilities using relative risks (RRs)
  grav_probs <- get_dist_grav(RRs)
  
  # Initialize variables for the multi_probs and prob_censor calculations
  multi_probs <- numeric(grav_cats * 2)
  prob_censor <- 0
  
  # Loop over each gravidity category to calculate probabilities for censoring and malaria status
  for (i in 1:grav_cats) {
    # Probability of not having malaria and not being censored
    multi_probs[2 * i - 1] <- (1 - malaria_prevalence[i]) * grav_probs[i] * (1 - prob_cens_with_params$prob_censor_no_malaria[i])
    prob_censor <- prob_censor + prob_cens_with_params$prob_censor_no_malaria[i] * (1 - malaria_prevalence[i]) * grav_probs[i]
    
    # Probability of having malaria and not being censored
    multi_probs[2 * i] <- malaria_prevalence[i] * grav_probs[i] * (1 - prob_cens_with_params$prob_censor_malaria[i])
    prob_censor <- prob_censor + prob_cens_with_params$prob_censor_malaria[i] * malaria_prevalence[i] * grav_probs[i]
  }
  
  # Normalize the multi_probs to ensure they sum to 1
  multi_probs <- multi_probs / sum(multi_probs)
  
  # Calculate the likelihood of observing the proportions by gravidity and malaria status
  malaria_grav_likelihood <- dmultinom(grav_mal_numbers, prob = multi_probs, log = TRUE)
  
  # Calculate the likelihood of observing the censored data
  censored_likelihood <- dbinom(censored_data$censored_number, censored_data$total_number, prob_censor, log = TRUE)
  
  # Return the total likelihood by combining GA likelihood, malaria gravidity likelihood, and censored likelihood
  return(GA_likelihood + malaria_grav_likelihood + censored_likelihood)
}

# Function to calculate the likelihood of the hemoglobin levels (HB_data)
process_block <- function(HB_data, site_FE, y_spline, G_non_infect, mal_spline, HB_sigma, grav_impact_df, cutoff = NULL) {
  # Calculate the mean hemoglobin levels based on various factors
  get_mean <- site_FE + y_spline[HB_data$adj_gestage] + G_non_infect[HB_data$gravidity] + mal_spline[HB_data$adj_gestage] * HB_data$malaria * grav_impact_df[HB_data$gravidity]
  
  # Calculate the likelihood using the normal distribution
  likelihood <- sum(dnorm(HB_data$hb_level, mean = get_mean, sd = HB_sigma, log = TRUE))
  
  # If a cutoff is provided, adjust the likelihood to account for censoring
  if (!is.null(cutoff)) {
    likelihood <- likelihood - sum(pnorm(cutoff, mean = get_mean, sd = HB_sigma, lower.tail = FALSE, log = TRUE))
  }
  
  # Return the likelihood
  return(likelihood)
}

# Function to calculate the log-likelihood for the censored dataset
r_loglike_w_censoring <- function(params, data, misc) { 
  country <- 1  # Define the country index (assuming only one country here)
  
  # Extract relevant parameters from the params vector
  HB_sigma <- params[grep("HB_sigma", names(params))]
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  G_non_infect <- c(0, params[grep("G[2-6]_non_infect", names(params))])
  y_knots <- params[grep("y_knot", names(params))]
  mal_knots <- params[grep("mal_knot", names(params))]
  malaria_prevalence <- plogis(params[grep("log_odds_malaria_prevalence_G[1-6]", names(params))])
  RR_Gs <- exp(params[grep("log_RR_G[2-6]", names(params))])
  shape_1_GA <- params[grep("shape_1_GA", names(params))]
  shape_2_GA <- params[grep("shape_2_GA", names(params))]
  
  # Define the gestational age range
  min <- misc$gestage_min
  max <- misc$gestage_max + 1
  
  # Generate the y_spline and mal_spline based on the cubic spline function
  y_spline <- cubic_spline(x = seq(min, max, length.out = 3), y = y_knots, x_pred = seq(min, max, length.out = (max - min) + 1))
  mal_spline <- cubic_spline(x = seq(min, max, length.out = 3), y = mal_knots, x_pred = seq(min, max, length.out = (max - min) + 1))
  
  # Calculate the immunity vector using the Hill function
  imm_vect <- hill_func(0:20, shape_hill, scale_hill)
  
  # Calculate the impact of gravidity on malaria effect
  grav_impact_df <- get_weight_impact_df(inf_history = misc$list_inf_histories[[country]], primi_prev = malaria_prevalence[1], imm_vect)
  
  # Initialize the total likelihood
  ret <- 0
  
  # Calculate the likelihood of the hemoglobin levels
  likelihood <- process_block(HB_data = data$country_df, site_FE = 0, y_spline = y_spline, G_non_infect = G_non_infect, mal_spline = mal_spline, grav_impact_df = grav_impact_df, HB_sigma = HB_sigma, cutoff = misc$cutoff)
  
  # Add the likelihood from the censoring process
  ret <- ret + process_censor_block(counts_df = data$counts_df, censored_data = data$censored_data, min_day = misc$gestage_min, max_day = max, cutoff = misc$cutoff, max_malaria = 1, site_FE = 0, y_spline = y_spline, mal_spline = mal_spline, HB_sigma = HB_sigma, G_non_infect = G_non_infect, grav_impact_df = grav_impact_df, grav_cats = max(data$counts_df$gravidity), shape_1_GA = shape_1_GA, shape_2_GA = shape_2_GA, malaria_prevalence = malaria_prevalence, RRs = RR_Gs, grav_mal_numbers = data$grav_mal_counts$total)
  
  # Add the hemoglobin level likelihood to the total likelihood
  ret <- ret + likelihood
  
  # Return the total likelihood
  return(ret)
}

# Function to calculate the log-likelihood for the non-censored dataset
r_loglike <- function(params, data, misc) { 
  country <- 1  # Define the country index (assuming only one country here)
  
  # Extract relevant parameters from the params vector
  HB_sigma <- params[grep("HB_sigma", names(params))]
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  G_non_infect <- c(0, params[grep("G[2-6]_non_infect", names(params))])
  y_knots <- params[grep("y_knot", names(params))]
  mal_knots <- params[grep("mal_knot", names(params))]
  PG_prev <- plogis(params[grep("log_odds_malaria_prevalence_G1", names(params))])
  
  # Define the gestational age range
  min <- misc$gestage_min
  max <- misc$gestage_max + 1
  
  # Generate the y_spline and mal_spline based on the cubic spline function
  y_spline <- cubic_spline(x = seq(min, max, length.out = 3), y = y_knots, x_pred = seq(min, max, length.out = (max - min) + 1))
  mal_spline <- cubic_spline(x = seq(min, max, length.out = 3), y = mal_knots, x_pred = seq(min, max, length.out = (max - min) + 1))
  
  # Calculate the immunity vector using the Hill function
  imm_vect <- hill_func(0:20, shape_hill, scale_hill)
  
  # Calculate the impact of gravidity on malaria effect
  grav_impact_df <- get_weight_impact_df(inf_history = misc$list_inf_histories[[country]], primi_prev = PG_prev, imm_vect)
  
  # Initialize the total likelihood
  ret <- 0
  
  # Calculate the likelihood of the hemoglobin levels
  likelihood <- process_block(HB_data = data$country_df, site_FE = 0, y_spline = y_spline, G_non_infect = G_non_infect, mal_spline = mal_spline, grav_impact_df = grav_impact_df, HB_sigma = HB_sigma, cutoff = NULL)
  
  # Calculate the likelihood of observing the prevalence in primigravid women
  prev_likelihood <- dbinom(data$prev_data$positive, data$prev_data$total, PG_prev, log = TRUE)
  
  # Add the hemoglobin level likelihood and prevalence likelihood to the total likelihood
  ret <- ret + likelihood + prev_likelihood
  
  # Return the total likelihood
  return(ret)
}

  