
# Function to compute the probability of censoring based on various factors
prob_cens <- function(counts_df,min_day, max_day, cutoff, max_malaria,site_FE, y_spline, mal_spline, HB_sigma, G_non_infect, grav_impact_df, grav_cats, shape_1_GA, shape_2_GA) {
  
  # Number of days
  n_days <- max_day - min_day + 1

  # Initialize the one-dimensional vector
  GA_beta_adj_censor <- numeric(n_days * (max_malaria + 1) * grav_cats)
  
  # Initialize the probabilities for being censored by infection status and gravidity
  prob_censor_no_malaria <- rep(0,grav_cats)
  prob_censor_malaria <- rep(0,grav_cats)
  base=G_non_infect[1:grav_cats]+site_FE
 # print(base)
  # Compute probability of enrollment for each day
  prob_GA <- pbeta(((seq(min_day, max_day) + 1 - min_day) / n_days), shape1 = shape_1_GA, shape2 = shape_2_GA) -
    pbeta((seq(min_day, max_day) - min_day) / n_days, shape1 = shape_1_GA, shape2 = shape_2_GA)
  #print(prob_GA)
  # Probability of being censored based on gestational age, gravidity, and malaria status
  no_malaria_probs_matrix <- t(outer(base, seq_len(n_days), function(base, i) {
    pnorm(cutoff, mean = base + y_spline[i], sd = HB_sigma)
  }))
 # mean_no_malaria <- t(outer(G_non_infect[1:grav_cats], seq_len(n_days), function(base, i) {
  #   base + y_spline[i]
  #}))
 
  malaria_probs_matrix <- mapply(function(base, effect) {
    sapply(seq_len(n_days), function(i) {
      #print(i)
      #print(base)
      #print(effect)
      #print(base + y_spline[i] + mal_spline[i] * effect)
      pnorm(cutoff, mean = base + y_spline[i] + mal_spline[i] * effect, sd = HB_sigma)
    })
  }, base, grav_impact_df[1:grav_cats])
  #malaria_mean_matrix <- mapply(function(base, effect) {
   # sapply(seq_len(n_days), function(i) {
      #print(i)
      #print(base)
      #print(effect)
      #print(base + y_spline[i] + mal_spline[i] * effect)
    #  mal_spline[i] * effect
    #})
  #}, G_non_infect[1:grav_cats], grav_impact_df[1:grav_cats])
  #print(malaria_probs_matrix)
  # Calculate the probabilities of being censored
  for (i in seq_len(n_days)) {
    prob_censor_no_malaria <- prob_censor_no_malaria + no_malaria_probs_matrix[i, ] * prob_GA[i]
    prob_censor_malaria <- prob_censor_malaria + malaria_probs_matrix[i, ] * prob_GA[i]
    
    for (j in 1:grav_cats) {
      # Calculate the index for the one-dimensional vector
      index_no_malaria <- (i - 1) * (max_malaria + 1) * grav_cats + j
      index_malaria <- index_no_malaria + grav_cats
      
      GA_beta_adj_censor[index_no_malaria] <- (1 - no_malaria_probs_matrix[i, j]) * prob_GA[i]
      GA_beta_adj_censor[index_malaria] <- (1 - malaria_probs_matrix[i, j]) * prob_GA[i]
    }
  }
  
  # Normalize GA_beta_adj_censor using vectorized operations
  norm_no_malaria <- 1 - prob_censor_no_malaria
  norm_malaria <- 1 - prob_censor_malaria
  
  # Create indices for each category
  no_malaria_indices <- rep(seq(1, grav_cats), each = n_days) + (0:(n_days-1)) * (max_malaria + 1) * grav_cats
  malaria_indices <- no_malaria_indices + grav_cats
  
  # Normalize in one step
  GA_beta_adj_censor[no_malaria_indices] <- GA_beta_adj_censor[no_malaria_indices] / rep(norm_no_malaria, each = n_days)
  GA_beta_adj_censor[malaria_indices] <- GA_beta_adj_censor[malaria_indices] / rep(norm_malaria, each = n_days)
  
 #index <- with(uncensored_dist_malawi, 
  #              (gestday - 1) + 
   #               (malaria) * n_days + 
    #              (grav_cat - 1) * n_days * (max_malaria + 1) + 1
  #)
  #index <- with(counts_df, (gestday - 1) * (max_malaria + 1) * grav_cats + malaria * grav_cats + grav_cat)
 # print(counts_df$count)
  #print(GA_beta_adj_censor[counts_df$index])
  GA_likelihood<- sum(counts_df$count * log(GA_beta_adj_censor[counts_df$index]))
  
 
 # print(prob_censor_no_malaria)
  #print(prob_censor_malaria)
  return(list(
    GA_likelihood= GA_likelihood,
    prob_censor_no_malaria = prob_censor_no_malaria,
    prob_censor_malaria = prob_censor_malaria
  ))
}
