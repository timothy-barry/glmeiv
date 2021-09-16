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
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 1
#' pi <- 0.2; n <- 1000; B <- 500
#' m_offset <- g_offset <- NULL
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
                               run_unknown_theta_precomputation = FALSE) {
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
  if (run_unknown_theta_precomputation) {
    for (i in seq(1L, B)) {
      print(i)
      fam <- MASS::negative.binomial(NA) %>% augment_family_object()
      m_precomp <- run_glmeiv_precomputation(y = data_list[[i]]$m,
                                             covariate_matrix = covariate_matrix,
                                             offset = m_offset,
                                             fam = fam)
      g_precomp <- run_glmeiv_precomputation(y = data_list[[i]]$g,
                                             covariate_matrix = covariate_matrix,
                                             offset = g_offset,
                                             fam = fam)
      attr(data_list[[i]], "m_precomp") <- m_precomp; attr(data_list[[i]], "g_precomp") <- g_precomp
    }
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


#' Process GLM-EIV results simulatr
#'
#' Processes the results outputted by a GLM-EIV method for use in simulatr.
#'
#' @param em_fit the fitted GLM-EIV model
#' @param s data frame of fitted coefficients and inferences
#' @param dat the data
#' @param save_membership_probs_mult integer giving multiple to save membership probabilities
#' @param time execution time (in s)
#'
#' @return a processed, long data frame
process_glmeiv_results_simulatr <- function(em_fit, s, time, dat, save_membership_probs_mult) {
  membership_prob_spread <- compute_mean_distance_from_half(em_fit$posterior_perturbation_probs)
  n_approx_1 <- sum(em_fit$posterior_perturbation_probs > 0.85)
  n_approx_0 <- sum(em_fit$posterior_perturbation_probs < 0.15)
  # output result
  meta_df <- tibble::tibble(parameter = "meta",
                          target = c("converged", "membership_probability_spread",
                                     "n_approx_0", "n_approx_1", "time"),
                          value = c(em_fit$converged, membership_prob_spread, n_approx_0, n_approx_1, time))
  out <- rbind(tidyr::pivot_longer(s, cols = -parameter, names_to = "target"), meta_df)
  # if i is a multiple of 250, save the posterior membership probabilities
  i <- attr(dat, "i")
  if (!is.null(i) && (i - 1 + save_membership_probs_mult) %% save_membership_probs_mult == 0) {
    out <- rbind(out, data.frame(parameter = "meta",
                               target = "membership_prob",
                               value = em_fit$posterior_perturbation_probs))
  }
 return(out)
}


#' Create simulatr specifier object (v2)
#'
#' @param param_grid grid of parameters
#' @param fixed_params list of fixed parameters
#' @param one_rep_times vector of times for each method, as well as the data generation procedure
#' @param methods a character vector giving the methods to run
#'
#' @return
#' @export
create_simulatr_specifier_object_v2 <- function(param_grid, fixed_params, one_rep_times, methods = c("glmeiv_slow", "glmeiv_fast", "thresholding")) {
  methods <- sort(methods)
  ####################################
  # 1. Define data_generator function.
  ####################################
  data_generator_object <- simulatr::simulatr_function(f = generate_full_data,
                                                       arg_names = c("m_fam", "m_intercept", "m_perturbation", "g_fam", "g_intercept", "g_perturbation", "pi",
                                                                     "n", "B", "covariate_matrix", "m_covariate_coefs", "g_covariate_coefs", "m_offset", "g_offset",
                                                                     "run_unknown_theta_precomputation"),
                                                       packages = "glmeiv",
                                                       loop = FALSE,
                                                       one_rep_time = one_rep_times[["generate_data_function"]])

  method_list <- c(
    if ("glmeiv_fast" %in% methods) simulatr::simulatr_function(f = run_glmeiv_at_scale_simulatr,
                                                                arg_names = c("m_fam", "g_fam", "covariate_matrix", "m_offset", "g_offset", "alpha",
                                                                              "n_em_rep", "save_membership_probs_mult", "pi_guess_range",
                                                                              "m_perturbation_guess_range", "g_perturbation_guess_range"),
                                                                packages = "glmeiv",
                                                                loop = TRUE,
                                                                one_rep_time = one_rep_times[["glmeiv_fast"]]) else NULL,
    if ("glmeiv_slow" %in% methods) simulatr::simulatr_function(f = run_glmeiv_random_init_simulatr,
                                                                arg_names = c("m_fam", "g_fam", "covariate_matrix", "m_offset",
                                                                              "g_offset", "alpha", "n_em_rep", "save_membership_probs_mult",
                                                                              "pi_guess_range", "m_perturbation_guess_range", "g_perturbation_guess_range",
                                                                              "m_intercept_guess_range", "g_intercept_guess_range", "m_covariate_coefs_guess_range",
                                                                              "g_covariate_coefs_guess_range"),
                                                                packages = "glmeiv",
                                                                loop = TRUE,
                                                                one_rep_time = one_rep_times[["glmeiv_slow"]]) else NULL,
    if ("thresholding" %in% methods) simulatr::simulatr_function(f = run_thresholding_method_simulatr,
                                                                 arg_names = c("g_intercept", "g_perturbation",
                                                                               "g_fam", "m_fam", "pi", "covariate_matrix",
                                                                               "g_covariate_coefs", "m_offset", "g_offset", "alpha"),
                                                                 packages = "glmeiv", loop = TRUE, one_rep_time = one_rep_times[["thresholding"]]) else NULL
    )
  names(method_list) <- methods

  ################################
  # 5. Instantiate sim_spec object
  ################################
  ret <- simulatr::simulatr_specifier(parameter_grid = param_grid,
                                      fixed_parameters = fixed_params,
                                      generate_data_function = data_generator_object,
                                      run_method_functions = method_list)
  return(ret)
}
