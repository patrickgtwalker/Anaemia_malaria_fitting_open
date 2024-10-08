# Load infection histories for the site
site_inf_histories <- readRDS("./site_inf_histories.RDS")
list_inf_histories <- list(
  site_inf_histories$site_inf_history[[1]]
)

get_dist_grav <- function(RRs) {
  # Calculate the distribution of gravidity categories based on relative risks (RRs)
  numerator <- c(1, RRs)  # Add baseline category (G1) with RR = 1
  denom <- sum(numerator)  # Sum of the numerators for normalization
  return(numerator / denom)  # Return the normalized distribution
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

# Parameters for the Beta distribution used to generate gestational age
GA_shape_1 <- 1
GA_shape_2 <- 1

# Log relative risks (log_RR) - decreasing sample size for categories 1-5, slight increase for category 6
log_RR <- c(-0.1, -0.2, -0.3, -0.4, 0.2)

# Calculate the relative distribution of sample sizes for each gravidity category using log_RR
grav_dist <- get_dist_grav(exp(log_RR))

# Generate sample sizes for each gravidity category using a multinomial distribution
sample_sizes <- as.vector(rmultinom(1, base_sample_size, grav_dist))

# Define log odds for prevalence (decreasing across gravidity categories)
Primigravid_prevalence <- 0.5
log_odds_base <- log(Primigravid_prevalence / (1 - Primigravid_prevalence))  # Calculate the base log odds
log_odds <- log_odds_base + c(0, -0.2, -0.4, -0.6, -0.8, -1.0)  # Adjust log odds for each gravidity category

# Convert log odds to prevalence
prevalence_by_grav <- exp(log_odds) / (1 + exp(log_odds))

# Define standard deviation for hemoglobin levels (HB_sigma)
HB_sigma <- 1.8

# Define knot points for y_spline and mal_spline
y_knots <- c(12, 10.5, 9.5)
mal_knots <- c(-0.5, -1.5, -1.25)

# Generate y_spline and mal_spline values based on the cubic spline function
y_spline <- cubic_spline(x = seq(gestage_min, gestage_max, length.out = 3), y = y_knots, x_pred = seq(gestage_min, gestage_max, length.out = (gestage_max - gestage_min) + 1))
mal_spline <- cubic_spline(x = seq(gestage_min, gestage_max, length.out = 3), y = mal_knots, x_pred = seq(gestage_min, gestage_max, length.out = (gestage_max - gestage_min) + 1))

# Create a data frame for the simulated malaria spline
simmed_mal_spline <- data.frame(
  gestage = seq(gestage_min, gestage_max, length.out = (gestage_max - gestage_min) + 1),
  mal_effect = -mal_spline  # Negate the malaria effect
)

# Hill function parameters for immunity vector
hill_shape <- 1
hill_scale <- 3

# Function to compute immunity vector using the Hill function
hill_func <- function(i, shape, scale) {
  return(1 / (1 + (i / scale)^shape))
}

# Generate immunity vector based on gravidity categories
imm_vect <- hill_func(0:20, hill_shape, hill_scale)
imm_effects <- get_weight_impact_df(site_inf_histories$site_inf_history[[1]], Primigravid_prevalence, imm_vect)

# Create a data frame for the simulated malaria immunity
simmed_mal_immunity <- data.frame(
  gravidity = 1:6,
  mal_effect = -mal_knots[2] * imm_effects  # Negate and scale malaria effects by immunity
)

# Non-infection effects (G_non_infect) for each gravidity category
G_non_infect_params <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

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
    country = 1,
    gravidity = i,
    gestage = gestage,
    adj_gestage = adj_gestage,
    malaria = malaria_status,
    hb_level = hb_levels
  )
}

# Combine all gravidity categories into one data frame
synthetic_data_df <- do.call(rbind, synthetic_data)

# Calculate data on overall risk in primigravidae for computational efficiency
PG_prev <- synthetic_data_df %>%
  filter(gravidity == 1) %>%
  summarise(positive = sum(malaria), total = n())

# Create a list of the different datasets for non-censored data
country_data_list_noncensor <- list(
  country_df = synthetic_data_df,
  prev_data = PG_prev
)

# Mimic censored data where women with HB < 7 were excluded
cutoff <- 7
synthetic_data_df_censor <- synthetic_data_df %>% filter(hb_level > cutoff)

# Provide data on the number of women censored
censored_number <- as.numeric(synthetic_data_df %>% filter(hb_level <= cutoff) %>% summarise(n()))
uncensored_number <- nrow(synthetic_data_df_censor)
censored_data <- list(
  censored_number = censored_number,
  total_number = censored_number + uncensored_number
)

# Pre-calculate the observed gestational age distribution in uncensored women
uncensored_dist <- synthetic_data_df_censor %>%
  group_by(adj_gestage, malaria, gravidity, .drop = FALSE) %>%
  summarise(count = n())

# Index based on gravidity and malaria status to align with function in prob_censor.R
uncensored_dist$index <- with(uncensored_dist, (adj_gestage - 1) * (max(uncensored_dist$malaria) + 1) * max(uncensored_dist$gravidity) + malaria * max(uncensored_dist$gravidity) + gravidity)

# Observed distribution by gravidity and malaria status
grav_mal_counts <- synthetic_data_df_censor %>% 
  group_by(gravidity, malaria) %>%
  summarise(total = n(), prop = n() / uncensored_number)

# Create a list of the different datasets for censored data
country_data_list_censor <- list(
  country_df = synthetic_data_df_censor,
  prev_data = PG_prev,
  grav_mal_counts = grav_mal_counts,
  counts_df = uncensored_dist,
  censored_data = censored_data
)

# Create a list of the different datasets for fitting without accounting for effect of censoring
country_data_list_wrong <- list(
  country_df = synthetic_data_df_censor,
  prev_data = PG_prev
)
