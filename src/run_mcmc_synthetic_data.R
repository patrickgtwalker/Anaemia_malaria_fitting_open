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
