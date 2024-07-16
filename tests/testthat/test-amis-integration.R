setwd("../..")

source("amis_integration.R")

sch_simulation <- get_amis_integration_package()
example_parameters <- sch_simulation$FixedParameters(
        # the higher the value of N, the more consistent the results will be
        # though the longer the simulation will take
        number_hosts = 10L,
        # no intervention
        coverage_file_name = "mansoni_coverage_scenario_0.xlsx",
        demography_name = "UgandaRural",
        # cset the survey type to Kato Katz with duplicate slide
        survey_type = "KK2",
        parameter_file_name = "mansoni_params.txt",
        coverage_text_file_storage_name = "Man_MDA_vacc.txt",
        # the following number dictates the number of events (e.g. worm deaths)
        # we allow to happen before updating other parts of the model
        # the higher this number the faster the simulation
        # (though there is a check so that there can't be too many events at once)
        # the higher the number the greater the potential for
        # errors in the model accruing.
        # 5 is a reasonable level of compromise for speed and errors, but using
        # a lower value such as 3 is also quite good
        min_multiplier = 5L
    )

test_that("Running the model should give us consistent results", {
    # Example prevalence map, with two locations, both with prevalence of 0.5
    prevalence_map <- matrix(c(0.5, 0.5), ncol = 1)

    tranmission_model <- build_transmission_model(prevalence_map, example_parameters, year_indices = c(23), 2)
    result <- tranmission_model(c(1L, 2L), matrix(c(3, 3, 0.3, 0.3), ncol = 2), 1)
    expect_equal(result, matrix(c(0.0, 0.5), ncol = 1), tolerance = 0.0)
})

test_that("Running the simulation on multiple time points gives multiple points back", {
    # Example prevalence map, with two locations, fitting to two time times
    # Both locations start at 0.031, and the second time point is 0.021
    prevalence_map <- vector("list", 2)
    prevalence_map[[1]]$data <- matrix(c(0.031, 0.031))
    prevalence_map[[2]]$data <- matrix(c(0.021, 0.021))

    year_indices <- c(0L, 23L)

    tranmission_model <- build_transmission_model(prevalence_map, example_parameters, year_indices, 2)
    result <- tranmission_model(c(1L, 2L), matrix(c(3, 3, 0.04, 0.04), ncol = 2), 1)
    expect_equal(result, matrix(c(0.1, 0.1, 0.1, 0.1), ncol = 2), tolerance = 0.5)
})

test_that("Running the simulation with different number of years specified compared to the prevalance map raises an error", {
    # Example prevalence map, with two locations, both with prevalence of 0.5
    prevalence_map <- matrix(c(0.5, 0.5), ncol = 1)

    year_indices <- c(0L, 23L)

    expect_error(build_transmission_model(prevalence_map, example_parameters, year_indices, 2), "Single time point prevalance map provided so should only request one year but 2 provided")
})

test_that("Running the simulation with different number of years specified compared to the prevalance map raises an error", {
    # Example prevalence map, with two locations, fitting to two time times
    # Both locations start at 0.031, and the second time point is 0.021
    prevalence_map <- vector("list", 2)
    prevalence_map[[1]]$data <- matrix(c(0.031, 0.031))
    prevalence_map[[2]]$data <- matrix(c(0.021, 0.021))

    year_indices <- c(23L)

    expect_error(build_transmission_model(prevalence_map, example_parameters, year_indices, 2), "Length of prevalance map \\(2\\) must match the number of years provided in year_indices \\(1\\)")
})

test_that("Running the AMIS integration on multiple time points should complete with the error about weight on particles", {
    # Example prevalence map, with two locations, fitting to two time times
    # Both locations start at 0.031, and the second time point is 0.021
    prevalence_map <- vector("list", 2)
    prevalence_map[[1]]$data <- matrix(c(0.031, 0.031))
    prevalence_map[[2]]$data <- matrix(c(0.021, 0.021))

    year_indices <- c(0L, 23L)

    #' the "dprior" function
    #' Note the second parameter _must_ be called log
    #' Unclear what this function does
    density_function <- function(parameters, log) {
        return(0.5)
    }

    #' The "rprior" function that returns a matrix whose columns are the parameters
    #' and each row is a sample
    rnd_function <- function(num_samples) {
        return(matrix(c(3, 0.04), ncol = 2, nrow = num_samples, byrow = TRUE))
    }

    prior <- list("dprior" = density_function, "rprior" = rnd_function)

    amis_params <- AMISforInfectiousDiseases::default_amis_params()
    amis_params$n_samples <- 2

    expect_error(AMISforInfectiousDiseases::amis(
        prevalence_map,
        build_transmission_model(prevalence_map, example_parameters, year_indices, 2),
        prior,
        amis_params
    ), "(No weight on any particles for locations in the active set.)|(the leading minor of order 2 is not positive definite)")
})

test_that("Running the AMIS integration should complete with the error about weight on particles", {
    # Example prevalence map, with two locations, both with prevalence of 0.5
    prevalence_map <- matrix(c(0.5, 0.5), ncol = 1)

    #' the "dprior" function
    #' Note the second parameter _must_ be called log
    #' Unclear what this function does
    density_function <- function(parameters, log) {
        return(0.5)
    }

    #' The "rprior" function that returns a matrix whose columns are the parameters
    #' and each row is a sample
    rnd_function <- function(num_samples) {
        return(matrix(c(3, 0.04), ncol = 2, nrow = num_samples, byrow = TRUE))
    }

    prior <- list("dprior" = density_function, "rprior" = rnd_function)

    amis_params <- AMISforInfectiousDiseases::default_amis_params()
    amis_params$n_samples <- 2

    expect_error(AMISforInfectiousDiseases::amis(
        prevalence_map,
        build_transmission_model(prevalence_map, example_parameters, year_indices = c(23), 2),
        prior,
        amis_params
    ), "(No weight on any particles for locations in the active set.)|(the leading minor of order 2 is not positive definite)")
})
