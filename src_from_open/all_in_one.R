"load_libraries.R"

library(easypackages)


libraries("foreign", "haven", "readstata13", "sas7bdat", "ggplot2", "dplyr", "tidyr", 
          "gridExtra", "ggnewscale", "ggpubr", "viridis", "drjacoby", 
          "matrixStats", "readxl", "grid","DescTools", "ggtext")

"cubic_spline.R"
# Function to compute cubic spline interpolation
cubic_spline <- function(x, y, x_pred) {
  # Number of intervals
  n <- length(x) - 1
  
  # Initialize arrays
  c <- l <- mu <- z <- rep(0, n + 1)
  h <- b <- d <- alpha <- rep(NA, n)
  
  # Compute the differences between consecutive x values
  for (i in 1:n) {
    h[i] <- x[i + 1] - x[i]
  }
  
  # Compute alpha values for the system of equations
  for (i in 2:n) {
    alpha[i] <- 3/h[i] * (y[i + 1] - y[i]) - 3/h[i - 1] * (y[i] - y[i - 1])
  }
  
  # Forward sweep to solve for c values
  l[1] <- 1
  for (i in 2:n) {
    l[i] <- 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
    mu[i] <- h[i]/l[i]
    z[i] <- (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
  }
  l[n + 1] <- 1
  
  # Backward substitution to calculate the coefficients
  for (i in n:1) {
    c[i] <- z[i] - mu[i] * c[i + 1]
    b[i] <- (y[i + 1] - y[i])/h[i] - h[i] * (c[i + 1] + 2 * c[i])/3
    d[i] <- (c[i + 1] - c[i])/(3 * h[i])
  }
  
  # Calculate the spline for each prediction point
  s <- rep(NA, length(x_pred))
  j <- 1
  for (i in seq_along(x_pred)) {
    while (x_pred[i] > x[j + 1]) {
      j <- j + 1
    }
    s[i] <- y[j] + b[j] * (x_pred[i] - x[j]) + c[j] * (x_pred[i] - x[j])^2 + d[j] * (x_pred[i] - x[j])^3
  }
  
  return(s)
}


"immune functions.R"
# Function to describe the decline in malaria effect using the Hill function
hill_func <- function(i, shape, scale) {
  return(1 / (1 + (i/scale)^shape))
}

# Function to compute weighted impact of gravidity category based on infection history
get_weight_impact_grav_cat <- function(inf_history_matrix, index, imm_vect_sub) {
  return(sum(imm_vect_sub * inf_history_matrix[index, ]))
}

# Optimized function to compute average impact on malaria effect from prior exposure across gravidity categories
get_weight_impact_df <- function(inf_history, primi_prev, imm_vect) {
  # Find the index that matches the primi_prev closest to inf_history
  index <- which.min(abs(inf_history$sim_details$primi_ANC_prev - primi_prev))
  # Calculate the weighted impacts for each gravidity category
  effects <- c(1,
               get_weight_impact_grav_cat(inf_history$inf_mat_g2, index, imm_vect[1:2]),
               get_weight_impact_grav_cat(inf_history$inf_mat_g3, index, imm_vect[1:3]),
               get_weight_impact_grav_cat(inf_history$inf_mat_g4, index, imm_vect[1:4]),
               get_weight_impact_grav_cat(inf_history$inf_mat_g5, index, imm_vect[1:5]),
               get_weight_impact_grav_cat(inf_history$inf_mat_g6, index, imm_vect[1:length(imm_vect)])
  )
  
  return(effects)
}

## prob_censor.R""
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


"generate_synthetic_data.R"

# Function to compute distribution of gravidity categories based on relative risks
get_dist_grav <- function(RRs) {
  numerator <- c(1, RRs)
  denom <- sum(numerator)
  return(numerator / denom)
}

# Set random seed for reproducibility
set.seed(56)

# Define the range for gestational age
gestage_min <- 80  # Minimum gestational age
gestage_max <- 200 # Maximum gestational age

# Number of pregnant women
base_sample_size <- 30000  # Total sample size across all gravidity categories

# Number of gravidity categories
grav_cats <- 6

GA_shape_1<-1
GA_shape_2<-1
# Log relative risks (log_RR) - decreasing sample size for categories 1-5, slight increase for category 6
log_RR <- c(-0.1, -0.2, -0.3, -0.4, 0.2)

# Calculate the relative distribution of sample sizes for each category using log_RR
grav_dist <- get_dist_grav(exp(log_RR))

# Generate sample sizes for each gravidity category using a multinomial distribution
sample_sizes <- as.vector(rmultinom(1, base_sample_size, grav_dist))

# Print the resulting sample sizes for each gravidity category
Primigravid_prevalence<-0.5

# Define log odds for prevalence (decreasing across gravidity categories)
log_odds_base <- log(Primigravid_prevalence / (1-Primigravid_prevalence))  # Base prevalence of 0.1
log_odds <- log_odds_base + c(0, -0.2, -0.4, -0.6, -0.8, -1.0)

prevalence_by_grav <- exp(log_odds) / (1 + exp(log_odds))


HB_sigma<-1.8
y_knots<-c(12,10.5,9.5)
mal_knots<-c(-0.5,-1.5,-1.25)



y_spline <- cubic_spline(x = seq(gestage_min,gestage_max,l=3), y = y_knots, x_pred = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1)))
mal_spline <- cubic_spline(x = seq(gestage_min,gestage_max,l=3), y = mal_knots, x_pred = seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1)))
simmed_mal_spline<-data.frame(gestage=seq(gestage_min,gestage_max,l=((gestage_max-gestage_min)+1)),
                              mal_effect=-mal_spline)

# Hill function parameters for immunity vector
hill_shape <-1
hill_scale <- 3

# Function to compute immunity vector using the Hill function
hill_func <- function(i, shape, scale) {
  return(1 / (1 + (i / scale)^shape))
}

# Immunity vector based on gravidity categories
imm_vect <- hill_func(0:20, hill_shape, hill_scale)
imm_effects<-get_weight_impact_df(site_inf_histories$site_inf_history[[1]],Primigravid_prevalence,imm_vect)

simmed_mal_immunity=data.frame(gravidity=1:6,
                               mal_effect=-mal_knots[2]*imm_effects)

# Non-infection effects (G_non_infect) for each gravidity category
G_non_infect_params <- c(0,0.1, 0.2, 0.3, 0.4, 0.5)


# Initialize empty lists to hold the synthetic data
synthetic_data <- list()

for (i in 1:grav_cats) {
  # Generate gestational ages for each individual
  adj_gestage<-floor(rbeta(n=sample_sizes[i],shape1 = GA_shape_1,shape2 = GA_shape_2)*(gestage_max-gestage_min))+1
  gestage <- adj_gestage+gestage_min-1
  
  
  # Generate malaria status (1 = malaria, 0 = no malaria) based on the prevalence
  malaria_status <- rbinom(sample_sizes[i], 1, prevalence_by_grav[i])
  
  # Access spline values based on rounded gestational age
  y_spline_values <- y_spline[adj_gestage]
  mal_spline_values <- mal_spline[adj_gestage]
  
  
  # Generate hemoglobin levels or other outcomes based on G_non_infect and malaria status
  hb_levels <- rnorm(sample_sizes[i], mean = G_non_infect_params[i] + y_spline_values + mal_spline_values *malaria_status*imm_effects[i], sd = HB_sigma)
  
  # Create a data frame for this gravidity category
  synthetic_data[[i]] <- data.frame(
    country=1,
    gravidity = i,
    gestage = gestage,
    adj_gestage=adj_gestage,
    malaria = malaria_status,
    hb_level = hb_levels
  )
}

# Combine all gravidity categories into one data frame
synthetic_data_df <- do.call(rbind, synthetic_data)
# precalculate data on overall risk in primigravidae for computational efficiency
PG_prev<-synthetic_data_df%>%filter(gravidity==1)%>%
  summarise(positive=sum(malaria),total=n())

## create a list of the different datasets 
country_data_list_noncensor<- list(country_df=synthetic_data_df,
                                   prev_data=PG_prev
)
### now lets censor the data to mimic studies where women with HB<7 where excluded
cutoff<-7
synthetic_data_df_censor<-synthetic_data_df%>%filter(hb_level>cutoff)

## provide data on # of women censored
censored_number<-as.numeric(synthetic_data_df%>%filter(hb_level<=cutoff)%>%summarise(n()))
uncensored_number<-nrow(synthetic_data_df_censor)
censored_data<-list(censored_number=censored_number,
                    total_number=censored_number+uncensored_number)

## pre-calculate the observed GA distribution in uncensored women for computational efficiency
uncensored_dist<-synthetic_data_df_censor%>%
  group_by(adj_gestage,malaria,gravidity,.drop=F)%>%
  summarise(count=n())

### index based on gravidity and malaria status to align with function in prob_censor.R for computational efficiency
uncensored_dist$index<-with(uncensored_dist, (adj_gestage - 1) * (max(uncensored_dist$malaria) + 1) * max(uncensored_dist$gravidity) + malaria * max(uncensored_dist$gravidity) + gravidity)

## observed distribution by gravidity and malaria status 
grav_mal_counts = synthetic_data_df_censor %>% dplyr::group_by(gravidity, malaria)%>%
  summarise(total=n(), prop=n()/uncensored_number)

## create a list of the different datasets 
country_data_list_censor <- list(country_df=synthetic_data_df_censor,
                                 prev_data=PG_prev,
                                 grav_mal_counts=grav_mal_counts,
                                 counts_df=uncensored_dist,
                                 censored_data=censored_data
)

## create a list of the different datasets for fitting without accounting for effect of censoring.
country_data_list_wrong <- list(country_df=synthetic_data_df_censor,
                                prev_data=PG_prev
)
"define_params_simulated_data.R"

### define a param_df for an uncensored dataset - don't need GA params and only Primigravid prevalence
max_GG<-6
df_params <- define_params(name = "HB_sigma", min = 0.1, max = 10,# removed GA_gradient
                           
                           name = "shape_hill", min = 0.01, max = 10, # Gravidity effect hill function
                           name = "scale_hill", min = 0.01, max = 10 # Gravidity effect hill function
                           #name = "shape_1_GA", min = 0.01, max = 10, 
                           #name = "shape_2_GA", min = 0.01, max = 25
)# made wider
y_knot_params <- define_params(name = sprintf("y_knot_%s", 1:3), min = 5, max = 15)
mal_knot_params <- define_params(name = sprintf("mal_knot_%s", 1:3), min = -10, max = 10)
PG_prev_param <- define_params(name="log_odds_malaria_prevalence_G1",  min = -4.595, max = 1.736)

df_params_non_infect<-define_params(name=sprintf("G%s_non_infect", 2:grav_cats),  min = -5, max = 5)
df_params_non_censor<-rbind(df_params,df_params_non_infect,y_knot_params,mal_knot_params,PG_prev_param)

### define a param_df for a censored dataset - additonal GA and prevalence params
prev_params <- define_params(name=sprintf("log_odds_malaria_prevalence_G%s", 1:grav_cats),  min = -4.595, max = 1.736)

RR_params <- define_params(name=sprintf("log_RR_G%s", 2:grav_cats),  min = -5, max = 5)

GA_params<-define_params(name = "shape_1_GA", min = 0.01, max = 10,
                         name = "shape_2_GA", min = 0.01, max = 10)


df_params_censor<-rbind(df_params,df_params_non_infect,y_knot_params,mal_knot_params,prev_params,RR_params,GA_params)


"log_likeli_synthetic_data.R"
process_censor_block <- function(counts_df,censored_data,
                                 min_day,max_day, cutoff, max_malaria,site_FE, y_spline, mal_spline, HB_sigma, G_non_infect,grav_impact_df,grav_cats,shape_1_GA, shape_2_GA,
                                 malaria_prevalence,RRs,grav_mal_numbers){
  
  prob_cens_with_params<-prob_cens(counts_df=counts_df,min_day=min_day,max_day=max_day, cutoff=cutoff, max_malaria=max_malaria, 
                                   site_FE=site_FE,y_spline=y_spline, mal_spline=mal_spline, HB_sigma=HB_sigma, G_non_infect=G_non_infect,grav_impact_df=grav_impact_df,
                                   grav_cats=grav_cats, shape_1_GA= shape_1_GA,shape_2_GA=shape_2_GA)
  
  
  
  ## likelihood of observed GA adjusting for likelihood it appears in the dataset given malaria status and gravidity
  # print(log(prob_cens_with_params$GA_beta_adj_censor[(floor(HB_data$gestage)-misc$gestage_min_wkenya+1)[1:2],(HB_data$malaria+1)[1:2],HB_data$grav_cat[1:2]]))    ### get the proportions by gravidity according to the RRs
  #print(log(grab_cens_GA_likeli(prob_cens_with_params$GA_beta_adj_censor,1,HB_data,misc$gestage_min_wkenya)))    ### get the proportions by gravidity according to the RRs
  #print(dim(prob_cens_with_params$GA_beta_adj_censor))
  #print(max(HB_data$newgestage))
  
  
  # GA_likelihood <-sum(log(sapply(1:length(HB_data$newgestage),grab_cens_GA_likeli,grav_likeli=prob_cens_with_params$GA_beta_adj_censor,data=HB_data,min_GA=min_day)))
  GA_likelihood<-prob_cens_with_params$GA_likelihood
  # Calculate gravidity probabilities
  grav_probs<-get_dist_grav(RRs)
  
  
  # Initialize variables
  multi_probs <- numeric(grav_cats*2)
  prob_censor <- 0
  
  
  # Loop over each gravidity category (G1 to G6)
  for (i in 1:grav_cats) {
    # Probability of not having malaria and not being censored
    multi_probs[2 * i - 1] <- (1 - malaria_prevalence[i]) * grav_probs[i] * (1 - prob_cens_with_params$prob_censor_no_malaria[i])
    prob_censor <- prob_censor + prob_cens_with_params$prob_censor_no_malaria[i] * (1 - malaria_prevalence[i]) * grav_probs[i]
    
    # Probability of having malaria and not being censored
    multi_probs[2 * i] <- malaria_prevalence[i] * grav_probs[i] * (1 - prob_cens_with_params$prob_censor_malaria[i])
    prob_censor <- prob_censor + prob_cens_with_params$prob_censor_malaria[i] * malaria_prevalence[i] * grav_probs[i]
  }
  
  # Normalize multi_probs to sum to 1
  multi_probs <- multi_probs / sum(multi_probs)
  ### get the likelihood of observed proportions by gravidity and malaria status
  
  malaria_grav_likelihood<-dmultinom(grav_mal_numbers, prob=multi_probs,log=TRUE)  #*`WKenya suffix here`
  censored_likelihood <- dbinom(censored_data$censored_number, censored_data$total_number, prob_censor, log=TRUE)
  
  return(GA_likelihood + malaria_grav_likelihood + censored_likelihood)
  
}


process_block <- function(HB_data, site_FE, y_spline, G_non_infect, mal_spline, HB_sigma,grav_impact_df, cutoff = NULL) {
  get_mean <- site_FE + y_spline[HB_data$adj_gestage] + G_non_infect[HB_data$gravidity]+ mal_spline[HB_data$adj_gestage]* HB_data$malaria*grav_impact_df[HB_data$gravidity]
  
  likelihood <- sum(dnorm(HB_data$hb_level, mean = get_mean, sd = HB_sigma, log = TRUE))
  
  if (!is.null(cutoff)) {
    likelihood <- likelihood - sum(pnorm(cutoff, mean = get_mean, sd = HB_sigma, lower.tail = FALSE, log = TRUE))
  }
  
  return(likelihood)
}


r_loglike_w_censoring <- function(params, data, misc) { 
  country<-1
  
  HB_sigma <- params[grep("HB_sigma", names(params))]
  
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  G_non_infect <- c(0, params[grep("G[2-6]_non_infect", names(params))])
  
  y_knots <- params[grep("y_knot", names(params))]
  mal_knots <- params[grep("mal_knot", names(params))]
  malaria_prevalence <- plogis(params[grep("log_odds_malaria_prevalence_G[1-6]", names(params))])
  RR_Gs <- exp(params[grep("log_RR_G[2-6]", names(params))])
  
  shape_1_GA <- params[grep("shape_1_GA", names(params))]
  shape_2_GA<- params[grep("shape_2_GA", names(params))]
  
  min = misc$gestage_min #50
  max = misc$gestage_max+1 #(238+1), bc want 238.5
  y_spline <- cubic_spline(x = seq(min,max,l=3), y = y_knots, x_pred = seq(min,max,l=((max-min)+1)))
  mal_spline <- cubic_spline(x = seq(min,max,l=3), y = mal_knots, x_pred = seq(min,max,l=((max-min)+1)))
  
  imm_vect<-hill_func(0:20,shape_hill,scale_hill)
  grav_impact_df <- get_weight_impact_df(inf_history = misc$list_inf_histories[[country]], #*`changed here`
                                         primi_prev = malaria_prevalence[1],  #*`changed here`
                                         imm_vect)
  
  ret = 0
  likelihood <- process_block(HB_data = data$country_df,
                              site_FE = 0,
                              y_spline = y_spline,
                              G_non_infect = G_non_infect,
                              mal_spline = mal_spline,
                              grav_impact_df = grav_impact_df,
                              HB_sigma = HB_sigma,
                              cutoff = misc$cutoff)
  
  
  # Prevalence and corresponding gravidity categories
  
  # this is the w kenya censoring block
  # wkenya censoring block is 9
  
  ret<-ret+process_censor_block(counts_df=data$counts_df,censored_data=data$censored_data,
                                min_day=misc$gestage_min,max_day=max, cutoff=misc$cutoff, max_malaria=1, 
                                site_FE=0,y_spline=y_spline, mal_spline=mal_spline, HB_sigma=HB_sigma, G_non_infect=G_non_infect,grav_impact_df=grav_impact_df,
                                grav_cats=max(data$counts_df$gravidity), shape_1_GA=shape_1_GA,shape_2_GA=shape_2_GA,
                                malaria_prevalence=malaria_prevalence,RRs=RR_Gs,grav_mal_numbers=data$grav_mal_counts$total)
  
  
  ret <- ret + likelihood
  return(ret)
}


r_loglike <- function(params, data, misc) { 
  country<-1
  
  HB_sigma <- params[grep("HB_sigma", names(params))]
  
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  G_non_infect <- c(0, params[grep("G[2-6]_non_infect", names(params))])
  
  y_knots <- params[grep("y_knot", names(params))]
  mal_knots <- params[grep("mal_knot", names(params))]
  
  PG_prev <- plogis(params[grep("log_odds_malaria_prevalence_G1", names(params))])
  
  
  min = misc$gestage_min #50
  max = misc$gestage_max+1 #(238+1), bc want 238.5
  y_spline <- cubic_spline(x = seq(min,max,l=3), y = y_knots, x_pred = seq(min,max,l=((max-min)+1)))
  mal_spline <- cubic_spline(x = seq(min,max,l=3), y = mal_knots, x_pred = seq(min,max,l=((max-min)+1)))
  
  imm_vect<-hill_func(0:20,shape_hill,scale_hill)
  grav_impact_df <- get_weight_impact_df(inf_history = misc$list_inf_histories[[country]], #*`changed here`
                                         primi_prev = PG_prev,  #*`changed here`
                                         imm_vect)
  ret = 0
  likelihood <- process_block(HB_data = data$country_df,
                              site_FE = 0,
                              y_spline = y_spline,
                              G_non_infect = G_non_infect,
                              mal_spline = mal_spline,
                              grav_impact_df = grav_impact_df,
                              HB_sigma = HB_sigma,
                              cutoff = NULL)
  
  prev_likelihood<-dbinom(data$prev_data$positive,data$prev_data$total,PG_prev,log=T)
  ret <- ret + likelihood+prev_likelihood
  return(ret)
}

"logprior.R"
r_logprior_noncensor <- function(params, misc) {
  HB_sigma <- params[grep("HB_sigma", names(params))]
  
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  
  
  return( dgamma(HB_sigma, shape = 0.01, rate = 0.01, log = TRUE) +
            dgamma(shape_hill,shape = 2, rate = 1, log = TRUE)+
            dgamma(scale_hill,shape = 2, rate = 1, log = TRUE)
  )
}

r_logprior_censor <- function(params, misc) {
  HB_sigma <- params[grep("HB_sigma", names(params))]
  
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  shape_1_GA <- params[grep("shape_1_GA", names(params))]
  shape_2_GA <- params[grep("shape_2_GA", names(params))]
  
  
  return( dgamma(HB_sigma, shape = 0.01, rate = 0.01, log = TRUE) +
            dgamma(shape_1_GA, shape = 0.01, rate = 0.01, log = TRUE) +
            dgamma(shape_2_GA, shape = 0.01, rate = 0.01, log = TRUE)+
            dgamma(shape_hill,shape = 2, rate = 1, log = TRUE)+
            dgamma(scale_hill,shape = 2, rate = 1, log = TRUE)
  )
}

"run_mcmc_synthetic_data.R"
site_inf_histories <- readRDS("./site_inf_histories.RDS")
list_inf_histories = list(
  site_inf_histories$site_inf_history[[1]]
)

source("./src/load_libraries.R")
source("./src/cubic_spline.R")
source("./src/immune_functions.R")
source("./src/prob_censor.R")
source("./src/generate_synthetic_data.R")
source("./src/define_params_simulated_data.R")
source("./src/log_likeli_synthetic_data.R")
source("./src/log_prior.R")


run_noncensored_data<-run_mcmc(data = country_data_list_noncensor,
                               df_params = df_params_non_censor,
                               misc = list(gestage_min = 80,
                                           gestage_max=200,
                                           cutoff = NULL,
                                           list_inf_histories = list_inf_histories
                               ),
                               loglike = r_loglike,
                               logprior = r_logprior_noncensor, # currently has no prior
                               burnin = 1000,
                               samples = 1000, 
                               pb_markdown = FALSE,
                               chains = 1
)

run_censored_data_correct<-run_mcmc(data = country_data_list_censor,
                                    df_params = df_params_censor,
                                    misc = list(gestage_min = 80,
                                                gestage_max=200,
                                                cutoff = 7,
                                                list_inf_histories = list_inf_histories
                                    ),
                                    loglike = r_loglike_w_censoring,
                                    logprior = r_logprior_censor, # currently has no prior
                                    burnin = 1000,
                                    samples = 1000, 
                                    pb_markdown = FALSE,
                                    chains = 1
)


run_censored_data_incorrect<-run_mcmc(data = country_data_list_wrong,
                                      df_params = df_params_non_censor,
                                      misc = list(gestage_min = 80,
                                                  gestage_max=200,
                                                  cutoff = NULL,
                                                  list_inf_histories = list_inf_histories
                                      ),
                                      loglike = r_loglike,
                                      logprior = r_logprior_noncensor, # currently has no prior
                                      burnin = 1000,
                                      samples = 1000, 
                                      pb_markdown = FALSE,
                                      chains = 1
)
