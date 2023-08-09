#' Generate full data
#'
#' Generates B copies of a full GLM-EIV dataset.
#'
#' @param m_fam family for m
#' @param m_intercept m intercept
#' @param m_perturbation m perturbation coefficient
#' @param g_fam family for g
#' @param g_intercept g intercept
#' @param g_perturbation g perturbation coefficient
#' @param pi probability of perturbation
#' @param n number of cells
#' @param B number of i.i.d. datasets to generate
#' @param covariate_matrix the (fixed) covariate matrix of confounding factors
#' @param m_covariate_coefs coefficients for technical factors in model for m
#' @param g_covariate_coefs coefficients for technical factors in model for g
#' @param m_offset optional (fixed) offset vector for m
#' @param g_offset optional (fixed) offset vector for g
#' @param run_unknown_theta_precomputation optional; if TRUE, runs the negative binomial (unknown theta) precomputation on all datasets and stores the results ("m_precomp," "g_precomp") as attributes.
#'
#' @return a list of length B of synthetic datasets with columns (p, m, g)
#' @export
#' @examples
#' \dontrun{
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' pi <- 0.2; n <- 1000; B <- 5
#' m_intercept <- log(0.01); m_perturbation <- log(0.25); g_intercept <- log(0.005); g_perturbation <- log(2.5)
#' m_offset <- log(rpois(n = n, lambda = 10000)); g_offset <- log(rpois(n = n, lambda = 5000))
#' # no covariates
#' m_covariate_coefs <- g_covariate_coefs <- covariate_matrix <- NULL
#' dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
#' g_intercept, g_perturbation, pi, n, B, covariate_matrix, m_covariate_coefs,
#' g_covariate_coefs, m_offset, g_offset)
#' # with covariates
#' covariate_matrix <- data.frame(p_mito = runif(n, 0, 10))
#' m_covariate_coefs <- -0.1; g_covariate_coefs <- 0.2
#' dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
#' g_intercept, g_perturbation, pi, n, B, covariate_matrix, m_covariate_coefs,
#' g_covariate_coefs, m_offset, g_offset)
#' }
generate_full_data <- function(m_fam, m_intercept, m_perturbation, g_fam, g_intercept, g_perturbation, pi, n,
                               B, covariate_matrix, m_covariate_coefs, g_covariate_coefs, m_offset, g_offset,
                               run_mrna_unknown_theta_precomputation = FALSE, run_grna_unknown_theta_precomputation = FALSE) {
  # if m_fam/g_fam in list form, extract
  if (!is(m_fam, "family")) m_fam <- m_fam[[1]]
  if (!is(g_fam, "family")) g_fam <- g_fam[[1]]

  # sample a B x n matrix of perturbation indicators
  perturbation_indicators <- matrix(data = stats::rbinom(n = n * B, size = 1, prob = pi), nrow = n, ncol = B)
  # call above for both m and g
  m_matrix <- generate_glm_data_sim(m_intercept, m_perturbation, perturbation_indicators, m_fam, covariate_matrix, m_covariate_coefs, m_offset, n, B)
  g_matrix <- generate_glm_data_sim(g_intercept, g_perturbation, perturbation_indicators, g_fam, covariate_matrix, g_covariate_coefs, g_offset, n, B)
  # Finally, create the data list
  data_list <- sapply(seq(1, B), function(i) {
    df <- data.frame(m = m_matrix[,i], g = g_matrix[,i], p = perturbation_indicators[,i])
    attr(df, "i") <- i
    return(df)
  }, simplify = FALSE)

  if (run_mrna_unknown_theta_precomputation) {
    if (m_fam$fam_str != "Negative Binomial") stop("Cannot run precomputation on mrna modality, as family is not NB.")
    for (i in seq(1L, B)) {
      fam <- MASS::negative.binomial(NA) |> augment_family_object()
      m_precomp <- run_glmeiv_precomputation(y = data_list[[i]]$m,
                                             covariate_matrix = covariate_matrix,
                                             offset = m_offset,
                                             fam = fam)
      attr(data_list[[i]], "m_precomp") <- m_precomp
    }
  }

  if (run_grna_unknown_theta_precomputation) {
    if (g_fam$fam_str != "Negative Binomial") stop("Cannot run precomputaton on gnra modality, as family is not NB.")
    fam <- MASS::negative.binomial(NA) |> augment_family_object()
    g_precomp <- run_glmeiv_precomputation(y = data_list[[i]]$g,
                                           covariate_matrix = covariate_matrix,
                                           offset = g_offset,
                                           fam = fam)
    attr(data_list[[i]], "g_precomp") <- g_precomp
  }

  return(data_list)
}


#' Generate glm data for simulation
#'
#' A lower-level function to generate GLM-EIV data.
#'
#' @param intercept intercept term
#' @param perturbation_coef coefficient for perturbation
#' @param perturbation_indicators an n times B matrix of perturbation indicators
#' @param fam family object describing response
#' @param covariate_matrix fixed matrix of covariates
#' @param covariate_coefs coefficients for technical factors
#' @param offset optional offset vector
#' @param n the number of examples
#' @param B the number of datasets to resample
#'
#' @export
#' @return an n times B matrix of synthetic response data
generate_glm_data_sim <- function(intercept, perturbation_coef, perturbation_indicators, fam, covariate_matrix, covariate_coefs, offset, n, B) {
  # compute theoretical conditional means
  conditional_means <- compute_theoretical_conditional_means(intercept, perturbation_coef, fam,
                                                             covariate_matrix, covariate_coefs, offset)
  mui0 <- conditional_means$mu0
  mui1 <- conditional_means$mu1
  varying_means <- length(mui0) >= 2
  # sample outputs
  y_matrix <- sapply(seq(1L, n), function(i) {
    row <- perturbation_indicators[i,]
    idx_mui0 <- which(row == 0)
    idx_mui1 <- which(row == 1)
    n_mui0 <- length(idx_mui0)
    n_mui1 <- length(idx_mui1)
    out <- numeric(length = B)
    out[idx_mui0] <- fam$simulate_n_times_given_mu(n = n_mui0, mu = if (varying_means) mui0[i] else mui0)
    out[idx_mui1] <- fam$simulate_n_times_given_mu(n = n_mui1, mu = if (varying_means) mui1[i] else mui1)
    return(out)
  })
  return(t(y_matrix))
}


#' Create simulatr specifier object (v2)
#'
#' @param param_grid grid of parameters
#' @param fixed_params list of fixed parameters
#' @param one_rep_times vector of times for each method, as well as the data generation procedure
#' @param methods a character vector giving the methods to run
#'
#' @return a simulatr specifier object for use in simulatr
#' @export
create_simulatr_specifier_object <- function(param_grid, fixed_params, methods = c("glmeiv_slow", "glmeiv_fast", "thresholding", "unimodal_mixture")) {
  methods <- sort(methods)
  ####################################
  # 0. Define the evaluation functions
  ####################################
  obtain_target <- function(output, target_name) {
    # first, check validity
    meta_block <- output |> dplyr::filter(parameter == "meta")
    if ("converged" %in% meta_block$target) {
      target_exists <- (meta_block |> dplyr::filter(target == "converged") |> dplyr::pull(value)) == 1
    } else if ("fit_attempted" %in%  meta_block$target) {
      target_exists <- (meta_block |> dplyr::filter(target == "fit_attempted") |> dplyr::pull(value)) == 1
    } else {
      stop("output not recognized")
    }
    if (target_exists) {
      if (target_name == "time") {
        ret <- meta_block |> dplyr::filter(target == "time") |> dplyr::pull("value")
      } else {
        ret <- output |> dplyr::filter(parameter == "m_perturbation" & target == target_name) |>
          dplyr::pull("value")
      }
    } else {
      ret <- NA
    }
    return(ret)
  }

  # bias
  compute_bias <- function(output, ground_truth) {
    estimate <- obtain_target(output, "estimate")
    if (is.na(estimate)) NA else (estimate - ground_truth)
  }

  # mse
  compute_mse <- function(output, ground_truth) {
    estimate <- obtain_target(output, "estimate")
    if (is.na(estimate)) NA else (estimate - ground_truth)^2
  }

  # coverage
  compute_coverage <- function(output, ground_truth) {
    confint_upper <- obtain_target(output, "confint_upper")
    confint_lower <- obtain_target(output, "confint_lower")
    if (is.na(confint_upper) || is.na(confint_lower)) {
      NA
    } else {
      ground_truth >= confint_lower && ground_truth <= confint_upper
    }
  }

  # ci width
  compute_ci_width <- function(output, ground_truth) {
    confint_upper <- obtain_target(output, "confint_upper")
    confint_lower <- obtain_target(output, "confint_lower")
    if (is.na(confint_upper) || is.na(confint_lower)) {
      NA
    } else {
      confint_upper - confint_lower
    }
  }

  # time
  compute_time <- function(output, ground_truth) {
    time <- obtain_target(output, "time")
    if (is.na(time)) NA else time
  }

  # rejection indicator
  compute_rejection_indicator <- function(output, ground_truth) {
    reject <- obtain_target(output, "p_value") < 0.05
    if (is.na(reject)) NA else reject
  }

  evaluation_functions <- list(bias = compute_bias, mse = compute_mse,
                               coverage = compute_coverage, ci_width = compute_ci_width,
                               time = compute_time, rejection_indicator = compute_rejection_indicator)

  ####################################
  # 1. Define data_generator function.
  ####################################
  data_generator_object <- simulatr::simulatr_function(f = generate_full_data,
                                                       arg_names = formalArgs(generate_full_data),
                                                       packages = "glmeiv",
                                                       loop = FALSE)

  method_list <- c(
    if ("glmeiv_fast" %in% methods) simulatr::simulatr_function(f = run_glmeiv_at_scale_simulatr,
                                                                arg_names = formalArgs(run_glmeiv_at_scale_simulatr)[-1],
                                                                packages = "glmeiv",
                                                                loop = TRUE),
    if ("glmeiv_slow" %in% methods) simulatr::simulatr_function(f = run_glmeiv_random_init_simulatr,
                                                                arg_names = formalArgs(run_glmeiv_random_init_simulatr)[-1],
                                                                packages = "glmeiv",
                                                                loop = TRUE),
    if ("thresholding" %in% methods) simulatr::simulatr_function(f = run_thresholding_method_simulatr,
                                                                 arg_names = formalArgs(run_thresholding_method_simulatr)[-1],
                                                                 packages = "glmeiv", loop = TRUE),
    if ("unimodal_mixture" %in% methods) simulatr::simulatr_function(f = run_umimodal_mixture_method_simulatr,
                                                                     arg_names = formalArgs(run_thresholding_method_simulatr)[-1],
                                                                     packages = "glmeiv", loop = TRUE)
    )
  names(method_list) <- methods

  ################################
  # 5. Instantiate sim_spec object
  ################################
  sim_spec <- simulatr::simulatr_specifier(parameter_grid = param_grid,
                                           fixed_parameters = fixed_params,
                                           generate_data_function = data_generator_object,
                                           run_method_functions = method_list,
                                           evaluation_functions = evaluation_functions)
  return(sim_spec)
}
