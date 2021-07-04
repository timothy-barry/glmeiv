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


#' Random initialization
#'
#' Output a vector of length n; n * pi of the entries are randomly set to 1, all others to 0.
#'
#' @param n length of output vector
#' @param pi fraction of entries to set to 1
#'
#' @return randomly initialized vector
#' @export
random_initialization <- function(n, pi) {
  out <- integer(n)
  out[sample(x = seq(1, n), size = floor(pi * n), replace = FALSE)] <- 1L
  return(out)
}


#' Run em algorithm for simulatr using optimal threshold initialization
#'
#' @param dat_list a list of data frames, each of which has columns m, g, p.
#' @param g_intercept the intercept for the gRNA model.
#' @param g_perturbation the perturbation coefficient for the gRNA model.
#' @param g_fam family object describing gRNA model.
#' @param m_fam family object describing mRNA model.
#' @param pi probability of perturbation
#' @param covariate_matrix data frame of technical factors; can be null
#' @param g_covariate_coefs technical factor coefficients for gRNA model
#' @param m_offset optional offset vector for mRNA model
#' @param g_offset optional offset vector for gRNA model
#' @param alpha (1-alpha)% CIs returned
#' @param n_em_rep number of EM algorithm runs to conduct
#' @param p_flip probability of flipping a given initialization weight.
#'
#' @return a data frame with columns parameter, target, value, and run_id
#' @export
#'
#' @examples
#' \dontrun{
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 3
#' pi <- 0.2; n <- 1000; B <- 5; alpha <- 0.95; n_em_rep <- 5; p_flip <- 0.01
#' m_offset <- g_offset <- NULL
#' m_covariate_coefs <- g_covariate_coefs <- covariate_matrix <- NULL
#' dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
#' g_intercept, g_perturbation, pi, n, B, covariate_matrix,
#' m_covariate_coefs, g_covariate_coefs, m_offset, g_offset)
#' run_em_algo_simulatr_optimal_thresh(dat_list, g_intercept, g_perturbation,
#' g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset,
#' alpha, n_em_rep, p_flip)
#' }
run_em_algo_simulatr_optimal_thresh <- function(dat_list, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha, n_em_rep, p_flip) {
  # first, obtain the optimal boundary
  bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
  n_datasets <- length(dat_list)
  n <- nrow(dat_list[[1]])
  res_list <- lapply(X = seq(1, n_datasets), FUN = function(i) {
    dat <- dat_list[[i]]
    g <- dat$g
    phat <- as.integer(g > bdy)
    if (all(phat == 1) || all(phat == 0)) { # just randomly initialize instead
      initial_Ti1_matrix <- replicate(n = n_em_rep, random_initialization(n, pi))
    } else {
      initial_Ti1_matrix <- cbind(phat, replicate(n_em_rep - 1, flip_weights(phat, p_flip)))
    }
    em_fit <- run_em_algo_multiple_inits(m = dat$m, g = dat$g, m_fam = m_fam, g_fam = g_fam,
                                         covariate_matrix = covariate_matrix, initial_Ti1_matrix = initial_Ti1_matrix,
                                         m_offset = m_offset, g_offset = g_offset, return_best = TRUE)
    s <- run_inference_on_em_fit(em_fit, alpha) %>% dplyr::rename("parameter" = "variable")
    tidyr::pivot_longer(s, cols = -parameter, names_to = "target") %>%
      dplyr::add_row(parameter = "information", target = "converged", value = em_fit$converged) %>%
      dplyr::mutate(run_id = i)
  })
  return(do.call(rbind, res_list))
}


#' Flip weights
#'
#' @param w a binary (0/1) vector
#' @param p the expected fraction of weights to flip
#'
#' @return a new binary vector with E(p) weights flipped.
#'
#' @examples
#' w <- rbinom(n = 100, size = 1, prob = 0.5)
#' p <- 0.1
#' glmeiv:::flip_weights(w, p)
flip_weights <- function(w, p) {
  out <- (w + stats::rbinom(n = length(w), size = 1, prob = p)) %% 2
  return(out)
}
