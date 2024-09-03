# Function to compute the probability of censoring based on various factors
prob_cens <- function(counts_df, min_day, max_day, cutoff, max_malaria, site_FE, y_spline, mal_spline, HB_sigma, G_non_infect, grav_impact_df, grav_cats, shape_1_GA, shape_2_GA) {
  
  n_days <- max_day - min_day + 1  # Calculate the number of days in the range
  GA_beta_adj_censor <- numeric(n_days * (max_malaria + 1) * grav_cats)  # Initialize a one-dimensional vector
  
  prob_censor_no_malaria <- rep(0, grav_cats)  # Initialize probabilities for no malaria
  prob_censor_malaria <- rep(0, grav_cats)     # Initialize probabilities for malaria
  base = G_non_infect[1:grav_cats] + site_FE  # Adjust the base effect by adding site-specific effects
  
  # Calculate the probability of enrollment for each gestational day
  prob_GA <- pbeta(((seq(min_day, max_day) + 1 - min_day) / n_days), shape1 = shape_1_GA, shape2 = shape_2_GA) -
    pbeta((seq(min_day, max_day) - min_day) / n_days, shape1 = shape_1_GA, shape2 = shape_2_GA)
  
  # Calculate the probability of being censored for uninfected individuals
  no_malaria_probs_matrix <- t(outer(base, seq_len(n_days), function(base, i) {
    pnorm(cutoff, mean = base + y_spline[i], sd = HB_sigma)
  }))
  
  # Calculate the probability of being censored for infected individuals
  malaria_probs_matrix <- mapply(function(base, effect) {
    sapply(seq_len(n_days), function(i) {
      pnorm(cutoff, mean = base + y_spline[i] + mal_spline[i] * effect, sd = HB_sigma)
    })
  }, base, grav_impact_df[1:grav_cats])
  
  # Calculate the adjusted probabilities of being censored for each gravidity category and malaria status
  for (i in seq_len(n_days)) {
    prob_censor_no_malaria <- prob_censor_no_malaria + no_malaria_probs_matrix[i, ] * prob_GA[i]
    prob_censor_malaria <- prob_censor_malaria + malaria_probs_matrix[i, ] * prob_GA[i]
    
    for (j in 1:grav_cats) {
      index_no_malaria <- (i - 1) * (max_malaria + 1) * grav_cats + j
      index_malaria <- index_no_malaria + grav_cats
      
      GA_beta_adj_censor[index_no_malaria] <- (1 - no_malaria_probs_matrix[i, j]) * prob_GA[i]
      GA_beta_adj_censor[index_malaria] <- (1 - malaria_probs_matrix[i, j]) * prob_GA[i]
    }
  }
  
  # Normalize the probabilities to ensure they sum to 1
  norm_no_malaria <- 1 - prob_censor_no_malaria
  norm_malaria <- 1 - prob_censor_malaria
  
  no_malaria_indices <- rep(seq(1, grav_cats), each = n_days) + (0:(n_days-1)) * (max_malaria + 1) * grav_cats
  malaria_indices <- no_malaria_indices + grav_cats
  
  GA_beta_adj_censor[no_malaria_indices] <- GA_beta_adj_censor[no_malaria_indices] / rep(norm_no_malaria, each = n_days)
  GA_beta_adj_censor[malaria_indices] <- GA_beta_adj_censor[malaria_indices] / rep(norm_malaria, each = n_days)
  
  # Compute the likelihood of the observed gestational age data
  GA_likelihood <- sum(counts_df$count * log(GA_beta_adj_censor[counts_df$index]))
  
  return(list(
    GA_likelihood = GA_likelihood,
    prob_censor_no_malaria = prob_censor_no_malaria,
    prob_censor_malaria = prob_censor_malaria
  ))
}
