
max_GG<-6
df_params <- define_params(name = "HB_sigma", min = 0.1, max = 10,# removed GA_gradient
                           
                             name = "shape_hill", min = 0.01, max = 10, # Gravidity effect hill function
                             name = "scale_hill", min = 0.01, max = 10 # Gravidity effect hill function
                             #name = "shape_1_GA", min = 0.01, max = 10, 
                             #name = "shape_2_GA", min = 0.01, max = 25
                           )# made wider
y_knot_params <- define_params(name = sprintf("y_knot_%s", 1:3), min = 5, max = 15)
mal_knot_params <- define_params(name = sprintf("mal_knot_%s", 1:3), min = -20, max = 20)
PG_prev_param <- define_params(name="log_odds_malaria_prevalence_G1",  min = -4.595, max = 1.736)

df_params_non_infect<-define_params(name=sprintf("G%s_non_infect", 2:grav_cats),  min = -5, max = 5)


df_params_non_censor<-rbind(df_params,df_params_non_infect,y_knot_params,mal_knot_params,PG_prev_param)

prev_params <- define_params(name=sprintf("log_odds_malaria_prevalence_G%s", 1:grav_cats),  min = -4.595, max = 1.736)

RR_params <- define_params(name=sprintf("log_RR_G%s", 2:grav_cats),  min = -10, max = 10)

GA_params<-define_params(name = "shape_1_GA", min = 0.01, max = 10,
                         name = "shape_2_GA", min = 0.01, max = 25)


df_params_censor<-rbind(df_params,df_params_non_infect,y_knot_params,mal_knot_params,prev_params,RR_params,GA_params)

