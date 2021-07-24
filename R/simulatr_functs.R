#' Create simulatr specifier object
#'
#' Creates a simulatr_specifier object to run a simulation.
#'
#' @param param_grid a grid of parameters giving the parameter settings
#' @param fixed_params a list of fixed parameters
#' @param covariate_sampler a list of functions to generate a covariate matrix
#' @param one_rep_times a named list giving the single rep time (either scalar or vector) of each method.
#'
#' @return a simulatr_specifier object
#' @export
create_simulatr_specifier_object <- function(param_grid, fixed_params, one_rep_times, covariate_sampler = NULL) {
  ############################################
  # 1. Create covariate_matrix (if necessary);
  # update the fixed_params list.
  ############################################
  if (!is.null(covariate_sampler)) {
    set.seed(1)
    covariate_matrix <- as.data.frame(lapply(covariate_sampler, function(f) f(fixed_params[["n"]])))
    fixed_params[["covariate_matrix"]] <- covariate_matrix
  }
  fixed_params[["m_fam"]] <- augment_family_object(fixed_params[["m_fam"]])
  fixed_params[["g_fam"]] <- augment_family_object(fixed_params[["g_fam"]])

  #######################################
  # 2. Define data_generator function and
  # corresponding data_generator simulatr
  # function object. Below, code some
  # checks of correctness in the comment.
  #######################################
  data_generator_object <- simulatr::simulatr_function(f = generate_full_data,
                                             arg_names = c("m_fam", "m_intercept", "m_perturbation", "g_fam", "g_intercept", "g_perturbation", "pi", "n",
                                                           "B", "covariate_matrix", "m_covariate_coefs", "g_covariate_coefs", "m_offset", "g_offset"),
                                             packages = "glmeiv",
                                             loop = FALSE,
                                             one_rep_time = one_rep_times[["generate_data_function"]])

  ######################################
  # 3. Define threshold estimator method
  ######################################
  thresholding_method_object <- simulatr::simulatr_function(f = run_thresholding_method_simulatr,
                                                            arg_names = c("g_intercept", "g_perturbation", "g_fam", "m_fam", "pi", "covariate_matrix",
                                                                          "g_covariate_coefs", "m_offset", "g_offset", "alpha"),
                                                            packages = "glmeiv",
                                                            loop = FALSE,
                                                            one_rep_time = one_rep_times[["thresholding"]])

  ###############################
  # 4. Define EM algorithm method
  ###############################
  em_method_object <- simulatr::simulatr_function(f = run_em_algo_mixture_init,
                                                  arg_names = c("g_fam", "m_fam", "covariate_matrix", "m_offset", "g_offset", "alpha", "n_em_rep", "sd", "save_membership_probs_mult", "lambda"),
                                                  packages = "glmeiv",
                                                  loop = TRUE,
                                                  one_rep_time = one_rep_times[["em"]])

  #############################
  # 5. Finally, instantiate the
  # simulatr_specifier object
  #############################
  ret <- simulatr::simulatr_specifier(parameter_grid = param_grid,
                               fixed_parameters = fixed_params,
                               generate_data_function = data_generator_object,
                               run_method_functions = list(thresholding = thresholding_method_object,
                                                           em = em_method_object))
  return(ret)
}


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
                               B, covariate_matrix, m_covariate_coefs, g_covariate_coefs, m_offset, g_offset) {
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
#' @return an n times B matrix of synthetic response data
generate_glm_data_sim <- function(intercept, perturbation_coef, perturbation_indicators, fam, covariate_matrix, covariate_coefs, offset, n, B) {
  # compute theoretical conditional means
  conditional_means <- compute_theoretical_conditional_means(intercept, perturbation_coef, fam,
                                                             covariate_matrix, covariate_coefs, offset)
  mui0 <- conditional_means$mu0
  mui1 <- conditional_means$mu1
  # sample outputs
  y_matrix <- sapply(seq(1, n), function(i) {
    row <- perturbation_indicators[i,]
    idx_mui0 <- which(row == 0)
    idx_mui1 <- which(row == 1)
    n_mui0 <- length(idx_mui0)
    n_mui1 <- length(idx_mui1)
    out <- numeric(length = B)
    out[idx_mui0] <- fam$simulate_n_times_given_mu(n = n_mui0, mu = if (is.null(covariate_matrix)) mui0 else mui0[i])
    out[idx_mui1] <- fam$simulate_n_times_given_mu(n = n_mui1, mu = if (is.null(covariate_matrix)) mui1 else mui1[i])
    return(out)
  })
  return(t(y_matrix))
}


#' Compute theoretical conditional means
#'
#' Computes the conditional means of the GLM-EIV model given the ground truth.
#'
#' @param intercept intercept term
#' @param perturbation_coef coefficient corresponding to perturbation
#' @param fam family object describing response distribution
#' @param covariate_matrix (optional) matrix of technical covariates
#' @param covariate_coefs (optional) coefficients corresponding to technical covariates
#' @param offset (optional) vector of offsets
#'
#' @return the scalar (or vector, if covariate_matrix is supplied) of conditional means
compute_theoretical_conditional_means <- function(intercept, perturbation_coef, fam, covariate_matrix = NULL, covariate_coefs = NULL, offset = NULL) {
  # augment family object and set offset to 0, if necessary
  if (is.null(offset)) offset <- 0
  # compute the (theoretical) conditional linear components
  if (is.null(covariate_matrix)) {
    li0 <- intercept + offset
    li1 <- li0 + perturbation_coef
  } else {
    form_str <- paste0("~", paste0(colnames(covariate_matrix), collapse = " + "), " + 0")
    m <- stats::model.matrix(stats::as.formula(form_str), covariate_matrix)
    li0 <- intercept + as.numeric((m %*% covariate_coefs)) + offset
    li1 <- li0 + perturbation_coef
  }
  # compute the (theoretical) conditional means
  mui0 <- fam$linkinv(li0)
  mui1 <- fam$linkinv(li1)
  return(list(mu0 = mui0, mu1 = mui1))
}


classify_estimates_em <- function(sim_res, spread_thresh = 0.1, approx_0_thresh = 60, approx_1_thresh = 60, g_pert_lower = -0.25, g_pert_upper = Inf, m_pert_lower = -Inf, m_pert_upper = 0.25, pi_lower = 0, pi_upper = 0.3) {
  out <- dplyr::filter(sim_res, method == "em") %>% dplyr::group_by(id) %>%
    dplyr::summarize(
      confident_output = (value[target == "converged"] == 1
                          & value[target == "membership_probability_spread"] > spread_thresh
                          & value[target == "n_approx_0"] >= approx_0_thresh
                          & value[target == "n_approx_1"] >= approx_1_thresh),
      g_perturbation = value[parameter == "g_perturbation" & target == "estimate"],
      m_perturbation = value[parameter == "m_perturbation" & target == "estimate"],
      pi = value[parameter == "pi" & target == "estimate"],
      plausible_estimates = (g_perturbation >= g_pert_lower &
                               g_perturbation <= g_pert_upper &
                               m_perturbation >= m_pert_lower &
                               m_perturbation <= m_pert_upper &
                               pi >= pi_lower & pi <= pi_upper),
      grid_row_id = grid_row_id[1]) %>%
    dplyr::mutate(classification = paste0(ifelse(confident_output, "confident", "unconfident"), "-",
                                          ifelse(plausible_estimates, "plausible", "implausible")) %>% factor(),
                  valid = confident_output & plausible_estimates)
  return(out)
}


#' Obtain valid IDs
#'
#' Obtains valid IDs given the output of a simulatr experiment.
#'
#' @param sim_res sim_res object
#' @param spread_thresh spread threshold
#' @param approx_0_thresh number cells approximately 0 threshold
#' @param approx_1_thresh number cells approximately 1 threshold
#' @param g_pert_lower minimum value for g_pert
#' @param g_pert_upper maximum value for g_pert
#' @param m_pert_lower minimum value for m_pert
#' @param m_pert_upper maximum value for m_pert
#' @param pi_lower minimum value for pi
#' @param pi_upper maximum value for pi
#'
#' @return a list of length two; (i) a data frame giving the classification of each EM run, and (ii) a character vector of valid IDs for both EM and thresholding
#' @export
obtain_valid_ids <- function(sim_res, spread_thresh = 0.1, approx_0_thresh = 50, approx_1_thresh = 50, g_pert_lower = -0.3, g_pert_upper = Inf, m_pert_lower = -Inf, m_pert_upper = 0.3, pi_lower = 0, pi_upper = 0.5) {
  em_df <- classify_estimates_em(sim_res, spread_thresh, approx_0_thresh, approx_1_thresh, g_pert_lower, g_pert_upper, m_pert_lower, m_pert_upper, pi_lower, pi_upper)
  valid_thresh_ids <- sim_res %>%
    dplyr::filter(method == "thresholding") %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(valid_idx = (value[target == "fit_attempted"] == 1)) %>%
    dplyr::filter(valid_idx) %>% dplyr::pull(id) %>% as.character()
  valid_em_ids <- dplyr::filter(em_df, valid) %>% dplyr::pull(id) %>% as.character()
  valid_ids <- c(valid_em_ids, valid_thresh_ids)
  return(list(em_classifications = em_df, valid_ids = valid_ids))
}
