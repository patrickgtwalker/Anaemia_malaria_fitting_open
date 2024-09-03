generate_df_params_non_censor <- function(
    grav_cats = 6, 
    HB_sigma_range = c(0.1, 10),
    hill_shape_range = c(0.01, 10),
    hill_scale_range = c(0.01, 10),
    y_knot_range = c(5, 15),
    mal_knot_range = c(-10, 10),
    PG_prev_range = c(-4.595, 1.736),
    non_infect_range = c(-5, 5)
) {
  df_params <- define_params(
    name = "HB_sigma", min = HB_sigma_range[1], max = HB_sigma_range[2],
    name = "shape_hill", min = hill_shape_range[1], max = hill_shape_range[2],  # Gravidity effect Hill function
    name = "scale_hill", min = hill_scale_range[1], max = hill_scale_range[2]   # Gravidity effect Hill function
  )
  
  # Define knot parameters for y_spline and mal_spline
  y_knot_params <- define_params(name = sprintf("y_knot_%s", 1:3), min = y_knot_range[1], max = y_knot_range[2])
  mal_knot_params <- define_params(name = sprintf("mal_knot_%s", 1:3), min = mal_knot_range[1], max = mal_knot_range[2])
  PG_prev_param <- define_params(name = "log_odds_malaria_prevalence_G1", min = PG_prev_range[1], max = PG_prev_range[2])
  
  # Define non-infection parameters for each gravidity category
  df_params_non_infect <- define_params(name = sprintf("G%s_non_infect", 2:grav_cats), min = non_infect_range[1], max = non_infect_range[2])
  
  # Combine all parameter sets for non-censored data
  df_params_non_censor <- rbind(df_params, df_params_non_infect, y_knot_params, mal_knot_params, PG_prev_param)
  
  return(df_params_non_censor)
}


generate_df_params_censor <- function(
    grav_cats = 6, 
    HB_sigma_range = c(0.2, 10),
    hill_shape_range = c(0.01, 10),
    hill_scale_range = c(0.01, 10),
    y_knot_range = c(5, 15),
    mal_knot_range = c(-5, 5),
    prev_range = c(-4.595, 1.736),
    RR_range = c(-5, 5),
    GA_shape_range = c(0.01, 10),
    non_infect_range = c(-5, 5)
) {
  df_params <- define_params(
    name = "HB_sigma", min = HB_sigma_range[1], max = HB_sigma_range[2],
    name = "shape_hill", min = hill_shape_range[1], max = hill_shape_range[2],  # Gravidity effect Hill function
    name = "scale_hill", min = hill_scale_range[1], max = hill_scale_range[2]   # Gravidity effect Hill function
  )
  
  # Define knot parameters for y_spline and mal_spline
  y_knot_params <- define_params(name = sprintf("y_knot_%s", 1:3), min = y_knot_range[1], max = y_knot_range[2])
  mal_knot_params <- define_params(name = sprintf("mal_knot_%s", 1:3), min = mal_knot_range[1], max = mal_knot_range[2])
  
  # Define prevalence parameters for each gravidity category
  prev_params <- define_params(name = sprintf("log_odds_malaria_prevalence_G%s", 1:grav_cats), min = prev_range[1], max = prev_range[2])
  
  # Define relative risk parameters for each gravidity category
  RR_params <- define_params(name = sprintf("log_RR_G%s", 2:grav_cats), min = RR_range[1], max = RR_range[2])
  
  # Define non-infection parameters for each gravidity category
  df_params_non_infect <- define_params(name = sprintf("G%s_non_infect", 2:grav_cats), min = non_infect_range[1], max = non_infect_range[2])
  
  # Define GA shape parameters
  GA_params <- define_params(
    name = "shape_1_GA", min = GA_shape_range[1], max = GA_shape_range[2],
    name = "shape_2_GA", min = GA_shape_range[1], max = GA_shape_range[2]
  )
  
  # Combine all parameter sets for censored data
  df_params_censor <- rbind(df_params, df_params_non_infect, y_knot_params, mal_knot_params, prev_params, RR_params, GA_params)
  
  return(df_params_censor)
}
