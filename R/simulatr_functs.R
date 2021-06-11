#' Create simulatr specifier object
#'
#' Creates a simulatr_specifier object to run a simulation.
#'
#' @param param_grid a grid of parameters giving the parameter settings
#' @param fixed_params a list of fixed parameters
#' @param covariate_sampler a list of functions to generate a covariate matrix
#'
#' @return a simulatr_specifier object
#' @export
#'
#' @examples
#' library(simulatr)
#' param_grid <- create_param_grid(varying_values = list(pi = seq(0.0, 0.5, 0.025),
#' m_perturbation = seq(0.0, -2, -0.1),
#' g_perturbation = seq(0.0, 2, 0.1)),
#' baseline_values = list(pi = 0.25, m_perturbation = -1, g_perturbation = 1))
#' fixed_params <- list(
#'  seed = 4,
#'  n = 1000,
#'  B = 2000,
#'  n_processors = 5,
#'  m_intercept = 2,
#'  g_intercept = 1,
#'  m_fam = poisson(),
#'  g_fam = poisson(),
#'  alpha = 0.95
#' )
#'
#' # simulation study 1: no covariates
#' covariate_sampler <- NULL
#'
#' # simulation study 2: covariates present
#' covariate_sampler <- list(lib_size = function(n) rpois(n, 10),
#'                           p_mito = function(n) runif(n, 0, 10))
#' fixed_params[["m_covariate_coefs"]] <- c(0.2, -0.1)
#' fixed_params[["g_covariate_coefs"]] <- c(0.1, -0.2)
create_simulatr_specifier_object <- function(param_grid, fixed_params, covariate_sampler = NULL) {
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
                                             loop = FALSE)
  # With covariates: dat_list <- generate_full_data(m_fam = fixed_params$m_fam, m_intercept = fixed_params$m_intercept, m_perturbation = 1, g_fam = fixed_params$g_fam, g_intercept = fixed_params$g_intercept, g_perturbation = 1, pi = 0.2, n = fixed_params$n, B = fixed_params$B, covariate_matrix = fixed_params$covariate_matrix, m_covariate_coefs = fixed_params$m_covariate_coefs, g_covariate_coefs = fixed_params$g_covariate_coefs, m_offset = NULL, g_offset = NULL)
  # m: fit <- glm(formula = m ~ p + lib_size + p_mito, family = fixed_params$m_fam, data = dplyr::mutate(fixed_params$covariate_matrix, dat[[1]]), offset = fixed_params$m_offset); coef(fit); c(fixed_params$m_intercept, 1, fixed_params$m_covariate_coefs)
  # g: fit <- glm(formula = g ~ p + lib_size + p_mito, family = fixed_params$g_fam, data = dplyr::mutate(fixed_params$covariate_matrix, dat[[1]]), offset = fixed_params$g_offset); coef(fit); c(fixed_params$g_intercept, 1, fixed_params$g_covariate_coefs)
  # Without covariates: dat <- generate_full_data(m_fam = fixed_params$m_fam, m_intercept = fixed_params$m_intercept, m_perturbation = 1, g_fam = fixed_params$g_fam, g_intercept = fixed_params$g_intercept, g_perturbation = 1, pi = 0.2, n = fixed_params$n, B = fixed_params$B, covariate_matrix = NULL, m_covariate_coefs = NULL, g_covariate_coefs = NULL, m_offset = NULL, g_offset = NULL)
  # fit <- glm(formula = m ~ p, family = fixed_params$m_fam, data = dat[[1]], offset = fixed_params$m_offset); coef(fit); c(fixed_params$m_intercept, 1)
  # fit <- glm(formula = g ~ p, family = fixed_params$g_fam, data = dat[[1]], offset = fixed_params$g_offset); coef(fit); c(fixed_params$g_intercept, 1)

  ######################################
  # 3. Define threshold estimator method
  ######################################
  thresholding_method_object <- simulatr::simulatr_function(f = run_thresholding_method_simulatr,
                                                            arg_names = c("g_intercept", "g_perturbation", "g_fam", "m_fam", "pi", "covariate_matrix",
                                                                          "g_covariate_coefs", "m_offset", "g_offset", "alpha"),
                                                            packages = "glmeiv",
                                                            loop = TRUE)
  #

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
generate_full_data <- function(m_fam, m_intercept, m_perturbation, g_fam, g_intercept, g_perturbation, pi, n,
                               B, covariate_matrix, m_covariate_coefs, g_covariate_coefs, m_offset, g_offset) {
  # sample a B x n matrix of perturbation indicators
  perturbation_indicators <- matrix(data = stats::rbinom(n = n * B, size = 1, prob = pi), nrow = n, ncol = B)
  # call above for both m and g
  m_matrix <- generate_glm_data_sim(m_intercept, m_perturbation, perturbation_indicators, m_fam, covariate_matrix, m_covariate_coefs, m_offset, n, B)
  g_matrix <- generate_glm_data_sim(g_intercept, g_perturbation, perturbation_indicators, g_fam, covariate_matrix, g_covariate_coefs, g_offset, n, B)
  # Finally, create the data list
  data_list <- sapply(seq(1, B), function(i) {
    data.frame(m = m_matrix[,i], g = g_matrix[,i], p = perturbation_indicators[,i])
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
