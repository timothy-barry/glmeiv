#' Generate data from model
#'
#' Generates data from the GLM-EIV model. The model treats the (unobserved) binary perturbations as random and all other
#' covariates as fixed.
#'
#' @param m_fam the augmented family object for M.
#' @param g_fam the augmented family object for G.
#' @param m_coef the coefficients of the GLM for M. The vector should be named; two of the names should be "intercept"
#' and "perturbation," and the other names should coincide with the column names of covariate_matrix.
#' @param g_coef the coefficients of the GLM for G. This vector should be named in a similar way to m_coef.
#' @param pi the probability of perturbation.
#' @param covariate_matrix the matrix of observed covariates; assumed to be fixed. Should be of class data.frame.
#' @param m_offset (optional) offsets for M.
#' @param g_offset (optional) offsets for G.
#' @param n (optional, required if covariate_matrix NULL) the number of samples to generate
#'
#' @return a list containing m, g, and p.
#'
#' @examples
#' \dontrun{
#' library(magrittr)
#' n <- 10000
#' m_fam <- augment_family_object(poisson())
#' g_fam <- augment_family_object(poisson())
#' m_offset <- 0
#' g_offset <- 0
#' pi <- 0.2
#' covariate_matrix <- data.frame(p_mito = runif(n = n, 0, 10),
#'                                lib_size = rpois(n = n, lambda = 1))
#' m_coef <- c(-1, -2, 1, 0.5)
#' g_coef <- c(-1, 2, 1, 0.5)
#' generated_data <- generate_data_from_model(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix)
#' # verify that we recover ground truth when p is known
#' covariate_matrix_full <- data.frame(perturbation = generated_data$p) %>%
#' dplyr::mutate(covariate_matrix)
#' fit_m <- glm(formula = generated_data$m ~ ., family = m_fam, data = covariate_matrix_full)
#' fit_g <- glm(formula = generated_data$g ~ ., family = g_fam, data = covariate_matrix_full)
#' }
generate_data_from_model <- function(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix, m_offset = NULL, g_offset = NULL, n = NULL) {
  # augment family objects, if necessary
  if (is.null(m_fam$augmented)) m_fam <- augment_family_object(m_fam)
  if (is.null(g_fam$augmented)) g_fam <- augment_family_object(g_fam)
  # set offsets to zero, if necessary
  if (is.null(m_offset)) m_offset <- 0; if (is.null(g_offset)) g_offset <- 0
  # verify column names ok
  check_col_names(covariate_matrix)
  # sample unobserved binary covariate p
  if (is.null(covariate_matrix)) {
    if (is.null(n)) stop("covariate_matrix is null. You must supply n.")
  } else {
    n <- nrow(covariate_matrix)
  }
  p <- stats::rbinom(n = n, size = 1, prob = pi)
  # append p to covariate matrix
  if (is.null(covariate_matrix)) {
    covariate_matrix_augmented <- data.frame(perturbation = p)
  } else {
    covariate_matrix_augmented <- dplyr::mutate(covariate_matrix, perturbation = p) %>% dplyr::select(perturbation, everything())
  }

  m <- simulate_glm_data(coefs = m_coef, fam = m_fam, offsets = m_offset, X = covariate_matrix_augmented)
  g <- simulate_glm_data(coefs = g_coef, fam = g_fam, offsets = g_offset, X = covariate_matrix_augmented)
  return(list(m = m, g = g, p = p))
}


#' Simulate GLM data
#'
#' @param coefs the named coefficient vector; one of the names should be "intercept;" the other names should match the column names of X.
#' @param fam the augmented family object
#' @param offsets the offset vector; can be 0.
#' @param X a data frame
#'
#' @return
#' data simulated from the GLM model
#' @noRd
simulate_glm_data <- function(coefs, fam, offsets, X) {
  formula <- paste0("~", paste0(colnames(X), collapse = " + ")) %>% stats::as.formula()
  cov_model <- stats::model.matrix(formula, data = X)

  # compute the linear portion of the model
  l_is <- as.numeric((cov_model %*% coefs)) + offsets
  mu_is <- fam$linkinv(l_is)
  y <- fam$simulate_from_mus(mu_is)
  return(y)
}


#' Get quick simulated data
#'
#' Quickly generates simulated data.
#'
#' @param n number of examples
#'
#' @return a list containing the entries m, g, p, pi, covariate_matrix, m_coef, and g_coef.
#' @export
get_quick_simulated_data <- function(n = 1000) {
  m_fam <- augment_family_object(stats::poisson())
  g_fam <- augment_family_object(stats::poisson())
  m_offset <- 0
  g_offset <- 0
  pi <- 0.2
  covariate_matrix <- data.frame(p_mito = stats::runif(n = n, 0, 10),
                                 lib_size = stats::rpois(n = n, lambda = 1))
  m_coef <- c(-2, 3, 0.75, -0.5)
  g_coef <- c(-2, 3, 1, 0.5)
  out <- generate_data_from_model(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix)
  out$pi <- pi
  out$covariate_matrix <- covariate_matrix
  out$m_coef <- m_coef
  out$g_coef <- g_coef
  out$m_fam <- m_fam
  out$g_fam <- g_fam
  return(out)
}


#' Run GLM-EIV with known p
#'
#' @param m m data vector
#' @param g g data vector
#' @param m_fam family object for m
#' @param g_fam family object for g
#' @param covariate_matrix the matrix of observed covariates
#' @param p the vector of p's; these typically are unobserved but will be used as the initial weights here
#' @param m_offset (default NULL) the vector of offsets for m
#' @param g_offset (default NULL) the vector of offsets for g
#' @param n_runs (default 10) number of EM algo runs
#' @param p_flip (default 0.15) expected fraction of initial weights to flip in each run
#' @param ep_tol (detault 0.1) tolerance threshold for EM convergence
#' @param max_it (default 50) maximum number of EM iterations (per run)
#' @param alpha (default 0.95) confidence interval level
#' @param reduced_output (default TRUE) return only the best EM run (as determined by log-likelihood)?
#'
#' @return the best EM run
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- get_quick_simulated_data(5000)
#' em_coefs <- run_glmeiv_known_p(m = dat$m, g = dat$g, m_fam = dat$m_fam, g_fam = dat$m_fam,
#' covariate_matrix = dat$covariate_matrix, p = dat$p, m_offset = NULL, g_offset = NULL)
#' # dat$m_coef; dat$g_coef; dat$pi
#' }
run_glmeiv_known_p <- function(m, g, m_fam, g_fam, covariate_matrix, p, m_offset = NULL, g_offset = NULL, n_runs = 5, p_flip = 0.15, ep_tol = 0.5 * 1e-4, max_it = 50, alpha = 0.95, reduced_output = TRUE) {
  initial_Ti1_matrix <- replicate(n_runs, expr = flip_weights(p, p_flip))
  em_runs <- run_em_algo_multiple_inits(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, ep_tol = ep_tol, max_it = max_it)
  if (reduced_output) {
    best_run <- select_best_em_run(em_runs)
    out <- run_inference_on_em_fit(best_run, alpha)
  } else {
    out <- em_runs
  }
  return(out)
}

# With covariates: dat <- generate_full_data(m_fam = fixed_params$m_fam, m_intercept = fixed_params$m_intercept, m_perturbation = 1, g_fam = fixed_params$g_fam, g_intercept = fixed_params$g_intercept, g_perturbation = 1, pi = 0.2, n = fixed_params$n, B = fixed_params$B, covariate_matrix = fixed_params$covariate_matrix, m_covariate_coefs = fixed_params$m_covariate_coefs, g_covariate_coefs = fixed_params$g_covariate_coefs, m_offset = NULL, g_offset = NULL)
# m: fit <- glm(formula = m ~ p + lib_size + p_mito, family = fixed_params$m_fam, data = dplyr::mutate(fixed_params$covariate_matrix, dat[[1]]), offset = fixed_params$m_offset); coef(fit); c(fixed_params$m_intercept, 1, fixed_params$m_covariate_coefs)
# g: fit <- glm(formula = g ~ p + lib_size + p_mito, family = fixed_params$g_fam, data = dplyr::mutate(fixed_params$covariate_matrix, dat[[1]]), offset = fixed_params$g_offset); coef(fit); c(fixed_params$g_intercept, 1, fixed_params$g_covariate_coefs)
# Without covariates: dat <- generate_full_data(m_fam = fixed_params$m_fam, m_intercept = fixed_params$m_intercept, m_perturbation = 1, g_fam = fixed_params$g_fam, g_intercept = fixed_params$g_intercept, g_perturbation = 1, pi = 0.2, n = fixed_params$n, B = fixed_params$B, covariate_matrix = NULL, m_covariate_coefs = NULL, g_covariate_coefs = NULL, m_offset = NULL, g_offset = NULL)
# fit <- glm(formula = m ~ p, family = fixed_params$m_fam, data = dat[[1]], offset = fixed_params$m_offset); coef(fit); c(fixed_params$m_intercept, 1)
# fit <- glm(formula = g ~ p, family = fixed_params$g_fam, data = dat[[1]], offset = fixed_params$g_offset); coef(fit); c(fixed_params$g_intercept, 1)
