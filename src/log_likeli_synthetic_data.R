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
  likelihood <- process_block(HB_data = country_data_list$country_df,
                              site_FE = 0,
                              y_spline = y_spline,
                              G_non_infect = G_non_infect,
                              mal_spline = mal_spline,
                              grav_impact_df = grav_impact_df,
                              HB_sigma = HB_sigma,
                              cutoff = NULL)
  
  prev_likelihood<-dbinom(country_data_list$prev_data$positive,country_data_list$prev_data$total,PG_prev,log=T)
  ret <- ret + likelihood+prev_likelihood
  return(ret)
}
  