get_amis_integration_package <- function() {
  importlib <- reticulate::import("importlib")
  sch_simulation <- reticulate::import("sch_simulation.amis_integration.amis_integration")
  importlib$reload(sch_simulation)
  return(sch_simulation)
}

build_transmission_model <- function(prevalence_map, fixed_parameters, year_indices, num_cores) {
  if(is.list(prevalence_map)) {
    if(length(prevalence_map) != length(year_indices)) {
      error_string <- sprintf("Length of prevalance map (%i) must match the number of years provided in year_indices (%i)", length(prevalence_map), length(year_indices))
      stop(error_string)
    }
  }
  else {
    if(length(year_indices) != 1) {
      error_string <- sprintf("Single time point prevalance map provided so should only request one year but %i provided", length(year_indices))
      stop(error_string)
    }
  }

  sch_simulation <- get_amis_integration_package()
  transmission_model <- function(seeds, params, n_tims) {
    output <- sch_simulation$run_model_with_parameters(
      # If year indices in just a single element, without as.array it will
      # automatically be converted into a scalar
      seeds, params, fixed_parameters, as.array(year_indices), as.integer(num_cores)
    )
    return(output)
  }

  return(transmission_model)
}
