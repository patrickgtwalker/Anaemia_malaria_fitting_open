  # Load necessary libraries and source functions
  source("./R/load_libraries.R")
  source("./R/cubic_spline.R")
  source("./R/immune_functions.R")
  source("./R/prob_censor.R")
  source("./R/data_generation_functions.R")
  source("./R/define_params_functions.R")
  source("./R/log_likeli.R")
  source("./R/log_prior.R")

list_inf_histories<-readRDS(here::here("inst", "extdata", "site_inf_histories.rds"))
selected_inf_history<-list_inf_histories$site_inf_history[[1]]
synthetic_data<-generate_synthetic_data(site_inf_history=selected_inf_history)


country_data_list_noncensor<-generate_country_data_list_noncensor(synthetic_data$synthetic_data_df)

df_params_non_censor<-generate_df_params_non_censor()


# Run MCMC for non-censored data
run_noncensored_data <- run_mcmc(
  data = country_data_list_noncensor,
  df_params = df_params_non_censor,
  misc = list(
    gestage_min = 80,
    gestage_max = 200,
    cutoff = NULL,
    list_inf_histories = list_inf_histories
  ),
  loglike = r_loglike,
  logprior = r_logprior_noncensor,
  burnin = 1000,
  samples = 1000,
  pb_markdown = FALSE,
  chains = 1
)

# Run MCMC for censored data with correct model
run_censored_data_correct <- run_mcmc(
  data = country_data_list_censor,
  df_params = df_params_censor,
  misc = list(
    gestage_min = 80,
    gestage_max = 200,
    cutoff = 7,
    list_inf_histories = list_inf_histories
  ),
  loglike = r_loglike_w_censoring,
  logprior = r_logprior_censor,
  burnin = 100,
  samples = 100,
  pb_markdown = FALSE,
  chains = 1
)

# Run MCMC for censored data with incorrect model (not accounting for censoring)
run_censored_data_incorrect <- run_mcmc(
  data = country_data_list_wrong,
  df_params = df_params_non_censor,
  misc = list(
    gestage_min = 80,
    gestage_max = 200,
    cutoff = NULL,
    list_inf_histories = list_inf_histories
  ),
  loglike = r_loglike,
  logprior = r_logprior_noncensor,
  burnin = 100,
  samples = 100,
  pb_markdown = FALSE,
  chains = 1
)

# save runs
#saveRDS(run_censored_data_incorrect, "./run_censored_incorrect.rds")
#saveRDS(run_censored_data_correct, "./run_censored_correct.rds")
#saveRDS(run_noncensored_data, "./run_noncensored_data.rds")
