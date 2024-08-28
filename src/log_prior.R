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
