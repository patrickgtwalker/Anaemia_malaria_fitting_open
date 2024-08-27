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
