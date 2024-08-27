source("./src/load_libraries.R")
#source("./src/load_data.R")
source("./src/cubic_spline.R")
source("./src/immune_functions.R")
source("./src/load_inf_histories.R")
source("./src/prob_censor.R")
source("./src/generate_synthetic_data.R")
source("./src/define_params_simulated_data.R")
source("./src/log_likeli_synthetic_data.R")
source("./src/log_prior.R")



list_inf_histories = list(
  site_inf_histories$site_inf_history[[1]]
)
quick_run<-run_mcmc(data = country_data_list,
                 df_params = df_params,
                 misc = list(gestage_min = 80,
                             gestage_max=200,
                             cutoff = 7,
                             list_inf_histories = list_inf_histories
                             ),
         loglike = r_loglike_w_censoring,
         logprior = r_logprior, # currently has no prior
         burnin = 1000,
         samples = 1000, 
         pb_markdown = FALSE,
         chains = 1
)

drjacoby::plot_density(quick_run)
head(country_data_list[[1]])
r_logprior <- function(params, misc) {
  return(0)
}
list_inf_histories[[1]]
short_run_low_censor<-run_mcmc(data = country_data_list,
         df_params = df_params_mal_site_FEs,
         misc = list(gestage_min_wkenya = min(wkenya_data$gestage),
                     gestage_max_wkenya = max(wkenya_data$gestage),
                     gestage_min_malawi = min(wkenya_data$gestage),
                     gestage_max_malawi = max(wkenya_data$gestage),#gestage max + 1 for range
                     cutoff = 7,
                     grav_mal_numbers_wkenya = grav_mal_numbers_wkenya,
                     grav_mal_numbers_malawi = grav_mal_numbers_malawi,
                     list_inf_histories = list_inf_histories),
         loglike = r_loglike,
         logprior = r_logprior2, # currently has no prior
         burnin = 1000,
         samples = 1000, 
         pb_markdown = FALSE,
         chains = 1) # takes 23 hours

short_run<-run_mcmc(data = country_data_list,
                    df_params = df_params_mal_site_FEs,
                    misc = list(gestage_min_wkenya = min(wkenya_data$gestage),
                                gestage_max_wkenya = max(wkenya_data$gestage),
                                gestage_min_malawi = min(wkenya_data$gestage),
                                gestage_max_malawi = max(wkenya_data$gestage),#gestage max + 1 for range
                                cutoff = 7,
                                grav_mal_numbers_wkenya = grav_mal_numbers_wkenya,
                                grav_mal_numbers_malawi = grav_mal_numbers_malawi,
                                list_inf_histories = list_inf_histories),
                    loglike = r_loglike,
                    logprior = r_logprior2, # currently has no prior
                    burnin = 1000,
                    samples = 1000, 
                    pb_markdown = FALSE,
                    chains = 1) 

saveRDS(short_run,"./sped_up_1000_1000.rds")
plot_density(short_run)

run_refact<-function(){
for(i in c(1:11)){
print(i)
print(r_loglike_refact(data = country_data_list,
                 params = df_params_mal_site_FEs,
                 misc = list(gestage_min_wkenya = min(wkenya_data$gestage),
                             gestage_max_wkenya = max(wkenya_data$gestage),
                             gestage_min_malawi = min(malawi_istp_data$gestage),
                             gestage_max_malawi = max(malawi_istp_data$gestage),#gestage max + 1 for range
                             cutoff = 7,
                             grav_mal_numbers_wkenya = grav_mal_numbers_wkenya,
                             grav_mal_numbers_malawi = grav_mal_numbers_malawi,
                             burkina_inf_history = burkina_inf_history,
                             gambia_inf_history = gambia_inf_history, 
                             ghana_inf_history = ghana_inf_history,
                             mali_inf_history = mali_inf_history,
                             kenya_inf_history = kenya_inf_history,
                             malawi_inf_history = malawi_inf_history,
                             tanzania_inf_history = tanzania_inf_history,
                             kenya_istp_inf_history = kenya_istp_inf_history,
                             malawi_istp_inf_history = malawi_istp_inf_history,block=i))
)
}
}
run_refact()
-12020.31

run_check<-function(){
  for(i in c(1:11)){
    
    print(
r_loglike_check(data = country_data_list,
                params = df_params_mal_site_FEs,
                misc = list(gestage_min_wkenya = min(wkenya_data$gestage),
                            gestage_max_wkenya = max(wkenya_data$gestage),
                            gestage_min_malawi = min(malawi_istp_data$gestage),
                            gestage_max_malawi = max(malawi_istp_data$gestage),#gestage max + 1 for range
                            cutoff = 7,
                            grav_mal_numbers_wkenya = grav_mal_numbers_wkenya,
                            grav_mal_numbers_malawi = grav_mal_numbers_malawi,
                            burkina_inf_history = burkina_inf_history,
                            gambia_inf_history = gambia_inf_history, 
                            ghana_inf_history = ghana_inf_history,
                            mali_inf_history = mali_inf_history,
                            kenya_inf_history = kenya_inf_history,
                            malawi_inf_history = malawi_inf_history,
                            tanzania_inf_history = tanzania_inf_history,
                            kenya_istp_inf_history = kenya_istp_inf_history,
                            malawi_istp_inf_history = malawi_istp_inf_history,block=i)))

}
}
run_check()
run_refact()

r_loglike_check_old(data = country_data_list,
                params = df_params_mal_site_FEs,
                misc = list(gestage_min_wkenya = min(wkenya_data$gestage),
                            gestage_max_wkenya = max(wkenya_data$gestage),
                            gestage_min_malawi = min(malawi_istp_data$gestage),
                            gestage_max_malawi = max(malawi_istp_data$gestage),#gestage max + 1 for range
                            cutoff = 7,
                            grav_mal_numbers_wkenya = grav_mal_numbers_wkenya,
                            grav_mal_numbers_malawi = grav_mal_numbers_malawi,
                            burkina_inf_history = burkina_inf_history,
                            gambia_inf_history = gambia_inf_history, 
                            ghana_inf_history = ghana_inf_history,
                            mali_inf_history = mali_inf_history,
                            kenya_inf_history = kenya_inf_history,
                            malawi_inf_history = malawi_inf_history,
                            tanzania_inf_history = tanzania_inf_history,
                            kenya_istp_inf_history = kenya_istp_inf_history,
                            malawi_istp_inf_history = malawi_istp_inf_history,block=10))


??microbenchmark
microbenchmark::microbenchmark(new=run_refact()
                               ,
                               old=run_check()                              
)
