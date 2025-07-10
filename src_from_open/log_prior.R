# Function to calculate the log-prior for non-censored data
r_logprior_noncensor <- function(params, misc) {
  # Extract relevant parameters from the params vector
  HB_sigma <- params[grep("HB_sigma", names(params))]
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  
  # Return the sum of the log-priors for each parameter
  return(
    dgamma(HB_sigma, shape = 0.01, rate = 0.01, log = TRUE) +  # Prior for HB_sigma
      dgamma(shape_hill, shape = 2, rate = 1, log = TRUE) +      # Prior for shape_hill
      dgamma(scale_hill, shape = 2, rate = 1, log = TRUE)        # Prior for scale_hill
  )
}

# Function to calculate the log-prior for censored data
r_logprior_censor <- function(params, misc) {
  # Extract relevant parameters from the params vector
  HB_sigma <- params[grep("HB_sigma", names(params))]
  shape_hill <- params[grep("shape_hill", names(params))]
  scale_hill <- params[grep("scale_hill", names(params))]
  shape_1_GA <- params[grep("shape_1_GA", names(params))]
  shape_2_GA <- params[grep("shape_2_GA", names(params))]
  
  # Return the sum of the log-priors for each parameter
  return(
    dgamma(HB_sigma, shape = 0.01, rate = 0.01, log = TRUE) +  # Prior for HB_sigma
      dgamma(shape_1_GA, shape = 0.01, rate = 0.01, log = TRUE) + # Prior for shape_1_GA
      dgamma(shape_2_GA, shape = 0.01, rate = 0.01, log = TRUE) + # Prior for shape_2_GA
      dgamma(shape_hill, shape = 2, rate = 1, log = TRUE) +       # Prior for shape_hill
      dgamma(scale_hill, shape = 2, rate = 1, log = TRUE)         # Prior for scale_hill
  )
}
