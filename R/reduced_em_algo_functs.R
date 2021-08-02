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
#' @param exp_m_offset exponentiated m offsets (i.e., fitted values from mRNA GLM)
#' @param exp_g_offset exponentiated g offsets (i.e., fitted values from gRNA GLM)
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
#' m_fam <- g_fam <- augment_family_object(poisson())
#' m_intercept <- 0
#' g_intercept <- 0
#' m_perturbation <- log(0.6)
#' g_perturbation <- log(1.4)
#' pi <- 0.01
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
#' exp_m_offset <- exp(m_offset)
#' exp_g_offset <- exp(g_offset)
#' fit_univariate <- run_reduced_em_algo(m, g, exp_m_offset, exp_g_offset,
#' m_pert_guess, g_pert_guess, pi_guess, m_fam, g_fam)
run_reduced_em_algo <- function(m, g, exp_m_offset, exp_g_offset, m_pert_guess, g_pert_guess, pi_guess, m_fam, g_fam, ep_tol = 0.5 * 1e-4, max_it = 50, min_it = 3) {
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
    e_step <- run_e_step_reduced(m = m, g = g, exp_m_offset = exp_m_offset, exp_g_offset = exp_g_offset,
                                 curr_m_perturbation = curr_m_perturbation, curr_g_perturbation = curr_g_perturbation,
                                 curr_pi = curr_pi, m_fam = m_fam, g_fam = g_fam)
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
    m_step <- run_m_step_univariate(m = m, g = g, exp_m_offset = exp_m_offset,
                                    exp_g_offset = exp_g_offset, curr_Ti1s = curr_Ti1s,
                                    m_fam = m_fam, g_fam = g_fam)
    curr_m_perturbation <- m_step$fit_m_perturbation
    curr_g_perturbation <- m_step$fit_g_perturbation
  }
  out <- list(m_perturbation = curr_m_perturbation, g_perturbation = curr_g_perturbation,
              pi = curr_pi, log_lik = curr_log_lik, converged = converged,
              n_iterations = iteration, log_liks = log_liks)
  return(out)
}


run_e_step_reduced <- function(m, g, exp_m_offset, exp_g_offset, curr_m_perturbation, curr_g_perturbation, curr_pi, m_fam, g_fam) {
  m_mus_pert1 <- exp(curr_m_perturbation) * exp_m_offset
  m_mus_pert0 <- exp_m_offset
  g_mus_pert1 <- exp(curr_g_perturbation) * exp_g_offset
  g_mus_pert0 <- exp_g_offset
  e_step <- run_e_step(m_fam = m_fam, g_fam = g_fam, m = m, g = g,
                       m_mus_pert0 = m_mus_pert0, m_mus_pert1 = m_mus_pert1,
                       g_mus_pert0 = g_mus_pert0, g_mus_pert1 = g_mus_pert1, fit_pi = curr_pi)
  return(e_step)
}


run_m_step_univariate <- function(m, g, exp_m_offset, exp_g_offset, curr_Ti1s, m_fam, g_fam) {
  # fit the models for m and g
  fit_m <- fit_model_univariate(m, curr_Ti1s, exp_m_offset, m_fam)
  fit_g <- fit_model_univariate(g, curr_Ti1s, exp_g_offset, g_fam)
  # return fitted parameters, as well as log-likelihood
  out <- list(fit_m_perturbation = fit_m, fit_g_perturbation = fit_g)
  return(out)
}


fit_model_univariate <- function(v, curr_Ti1s, exp_offset, fam) {
  p1 <- sum(curr_Ti1s * v)
  p2 <- sum(curr_Ti1s * exp_offset)
  fit_param <- log(p1) - log(p2)
  return(fit_param)
}
