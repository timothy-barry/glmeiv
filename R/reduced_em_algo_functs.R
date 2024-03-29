#' Run univariate poisson EM algo
#'
#' This function aims to estimate three parameters: m_perturbation, g_perturbation, and pi.
#'
#' The function takes offsets for the mRNA and gRNA modalities. It is assumed that all auxilliary information
#' (baseline expression level, library size, batch effect and other technical factors, etc.) has been distilled into
#' the offset vectors.
#'
#' The function requires initial estimates for the parameters m_perturbation, g_perturbation, and pi.
#'
#' @param m mRNA counts
#' @param g gRNA counts
#' @param m_fitted fitted values on linear scale of mRNA model
#' @param g_fitted fitted values on linear scale of gRNA model
#' @param m_pert_guess initial guess for m perturbation parameter
#' @param g_pert_guess initial guess for g perturbation parameter
#' @param pi_guess initial guess for pi
#' @param ep_tol relative tolerance limit for EM algorithm
#' @param max_it maximum number of iterations before giving up
#' @param m_fam family object describing mRNA distribution
#' @param g_fam family object describing gRNA distribution
#' @param min_it mininum number of iterations before declaring convergence
#'
#' @return a list containing (i) the vector of estimates, (ii) log-likelihood of the fitted model, and (iii) the posterior membership probabilities.
#' @export
#'
#' @examples
#' set.seed(4)
#' # NB response
#' m_fam <- g_fam <- augment_family_object(MASS::negative.binomial(20))
#' m_intercept <- 0
#' g_intercept <- 0
#' m_perturbation <- log(0.6)
#' g_perturbation <- log(1.4)
#' pi <- 0.02
#' n <- 100000
#' m_offset <- log(stats::rpois(n, 100))
#' g_offset <- log(stats::rpois(n, 50))
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept, m_perturbation = m_perturbation,
#' g_fam = g_fam, g_intercept = g_intercept, g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = NULL, m_covariate_coefs = NULL, g_covariate_coefs = NULL, m_offset = m_offset,
#' g_offset = g_offset)[[1]]
#' m_fit <- glm(formula = m ~ p + 0, family = m_fam, data = dat, offset = m_offset)
#' g_fit <- glm(formula = g ~ p + 0, family = g_fam, data = dat, offset = g_offset)
#' m_pert_guess <- log(runif(1, 0.1, 1.5))
#' g_pert_guess <- log(runif(1, 0.5, 10))
#' pi_guess <- runif(1, 0, 0.05)
#' m <- dat$m
#' g <- dat$g
#' m_fitted <- m_offset
#' g_fitted <- g_offset
#' fit_univariate <- run_reduced_em_algo(m, g, m_fitted, g_fitted,
#' m_pert_guess, g_pert_guess, pi_guess, m_fam, g_fam)
#'
#' # Gaussian
#' m_fam <- g_fam <- augment_family_object(gaussian())
#' m_intercept <- 0
#' g_intercept <- 0
#' m_perturbation <- 2
#' g_perturbation <- -1
#' pi <- 0.02
#' n <- 100000
#' m_offset <- log(stats::rpois(n, 100))
#' g_offset <- log(stats::rpois(n, 50))
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept, m_perturbation = m_perturbation,
#' g_fam = g_fam, g_intercept = g_intercept, g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = NULL, m_covariate_coefs = NULL, g_covariate_coefs = NULL, m_offset = m_offset,
#' g_offset = g_offset)[[1]]
#' m_fit <- glm(formula = m ~ p + 0, family = m_fam, data = dat, offset = m_offset)
#' g_fit <- glm(formula = g ~ p + 0, family = g_fam, data = dat, offset = g_offset)
#' m_pert_guess <- 3
#' g_pert_guess <- -2
#' pi_guess <- runif(1, 0, 0.05)
#' m <- dat$m
#' g <- dat$g
#' m_fitted <- m_offset
#' g_fitted <- g_offset
#' fit_univariate <- run_reduced_em_algo(m, g, m_fitted, g_fitted,
#' m_pert_guess, g_pert_guess, pi_guess, m_fam, g_fam)
run_reduced_em_algo <- function(m, g, m_fitted, g_fitted, m_pert_guess, g_pert_guess, pi_guess, m_fam, g_fam, ep_tol = 0.5 * 1e-4, max_it = 50, min_it = 3, use_mrna_modality = TRUE) {
  set.seed(4)
  # set some basic variables
  converged <- FALSE
  prev_log_lik <- -Inf
  n <- length(m)
  curr_m_perturbation <- m_pert_guess; curr_g_perturbation <- g_pert_guess; curr_pi <- pi_guess
  iteration <- 1L
  log_liks <- numeric()
  # Iterate through E and M steps until convergence.
  while (TRUE) {
    # E step
    e_step <- run_e_step_reduced(m = m, g = g, m_fitted = m_fitted, g_fitted = g_fitted,
                                 curr_m_perturbation = curr_m_perturbation, curr_g_perturbation = curr_g_perturbation,
                                 curr_pi = curr_pi, m_fam = m_fam, g_fam = g_fam, use_mrna_modality = use_mrna_modality)
    curr_Ti1s <- e_step$Ti1s
    if (any(is.na(curr_Ti1s)) || all(curr_Ti1s == 0)) {
      curr_log_lik <- -Inf
      break()
    }
    curr_log_lik <- e_step$log_lik
    curr_pi <- e_step$new_pi
    log_liks <- c(log_liks, curr_log_lik)

    # Check for convergence
    tol <- compute_tolerance(curr_log_lik, prev_log_lik)
    if (tol < ep_tol && iteration >= min_it) {
      converged <- TRUE; break()
    } else {
      prev_log_lik <- curr_log_lik
      iteration <- iteration + 1L
      if (iteration >= max_it) break()
    }

    # M step
    m_step <- run_m_step_univariate(m = m, g = g, m_fitted = m_fitted,
                                    g_fitted = g_fitted, curr_Ti1s = curr_Ti1s,
                                    m_fam = m_fam, g_fam = g_fam, use_mrna_modality = use_mrna_modality)
    curr_m_perturbation <- m_step$fit_m_perturbation
    curr_g_perturbation <- m_step$fit_g_perturbation
  }
  out <- list(m_perturbation = curr_m_perturbation, g_perturbation = curr_g_perturbation,
              pi = curr_pi, log_lik = curr_log_lik, converged = converged,
              n_iterations = iteration, log_liks = log_liks, curr_Ti1s = curr_Ti1s)
  return(out)
}


run_e_step_reduced <- function(m, g, m_fitted, g_fitted, curr_m_perturbation, curr_g_perturbation, curr_pi, m_fam, g_fam, use_mrna_modality = TRUE) {
  if (use_mrna_modality) {
    m_mus_pert0 <- m_fam$linkinv(m_fitted)
    m_mus_pert1 <- m_fam$linkinv(m_fitted + curr_m_perturbation)
  } else {
    m_mus_pert0 <- m_mus_pert1 <- NULL
  }

  g_mus_pert0 <- g_fam$linkinv(g_fitted)
  g_mus_pert1 <- g_fam$linkinv(g_fitted + curr_g_perturbation)

  e_step <- run_e_step(m_fam = m_fam, g_fam = g_fam, m = m, g = g,
                       m_mus_pert0 = m_mus_pert0, m_mus_pert1 = m_mus_pert1,
                       g_mus_pert0 = g_mus_pert0, g_mus_pert1 = g_mus_pert1,
                       fit_pi = curr_pi, use_mrna_modality = use_mrna_modality)
  return(e_step)
}


run_m_step_univariate <- function(m, g, m_fitted, g_fitted, curr_Ti1s, m_fam, g_fam, use_mrna_modality = TRUE) {
  # fit the models for m and g
  if (use_mrna_modality) {
    fit_m <- fit_model_univariate(m, curr_Ti1s, m_fitted, m_fam)
  } else {
    fit_m <- NULL
  }
  fit_g <- fit_model_univariate(g, curr_Ti1s, g_fitted, g_fam)
  # return fitted parameters, as well as log-likelihood
  out <- list(fit_m_perturbation = fit_m, fit_g_perturbation = fit_g)
  return(out)
}


fit_model_univariate <- function(v, curr_Ti1s, fitted, fam) {
  if (fam$fam_str %in% c("poisson", "Negative Binomial")) {
    p1 <- sum(curr_Ti1s * v)
    p2 <- sum(curr_Ti1s * exp(fitted))
    fit_param <- log(p1) - log(p2)
  } else {
    fit_param <- sum(curr_Ti1s * (v - fitted))/sum(curr_Ti1s)
  }
  return(fit_param)
}
