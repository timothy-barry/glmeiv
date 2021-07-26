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
#' @return a list containing (i) the vector of estimates, (ii) log-likelihood of the fitted model, and (iii) the posterior membership probabilities.
#' @export
#'
#' @examples
#' m_fam <- poisson() %>% augment_family_object()
#' g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 0
#' g_intercept <- 0
#' m_perturbation <- log(0.5) # -0.69
#' g_perturbation <- log(1.2) # 0.18
#' pi <- 0.01
#' n <- 100000
#' m_offset <- log(stats::rpois(n, 100))
#' g_offset <- log(stats::rpois(n, 50))
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept, m_perturbation = m_perturbation, g_fam = g_fam,
#'                           g_intercept = g_intercept, g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#'                           covariate_matrix = NULL, m_covariate_coefs = NULL, g_covariate_coefs = NULL, m_offset = m_offset, g_offset = g_offset)[[1]]
#' m_fit <- glm(formula = m ~ p + 0, family = m_fam, data = dat, offset = m_offset)
#' g_fit <- glm(formula = g ~ p + 0, family = g_fam, data = dat, offset = g_offset)
#' m_pert_guess <- log(runif(1, 0.1, 1.5))
#' g_pert_guess <- log(runif(1, 0.5, 10))
#' pi_guess <- runif(1, 0, 0.05)
#' m <- dat$m
#' g <- dat$g
#' exp_m_offset <- exp(m_offset)
#' exp_g_offset <- exp(g_offset)
#' fit_univariate <- run_univariate_poisson_em_algo(m, g, exp_m_offset, exp_g_offset, m_pert_guess, g_pert_guess, pi_guess)
run_univariate_poisson_em_algo <- function(m, g, exp_m_offset, exp_g_offset, m_pert_guess, g_pert_guess, pi_guess, ep_tol = 0.5 * 1e-4, max_it = 50) {
  # set some basic variables
  converged <- FALSE
  prev_log_lik <- -Inf
  n <- length(m)
  m_fam <- g_fam <- augment_family_object(poisson())
  curr_m_perturbation <- m_pert_guess; curr_g_perturbation <- g_pert_guess; curr_pi <- pi_guess
  iteration <- 0L

  # iterate through E and M steps until convergence
  while (!converged) {
    # E step
    curr_Ti1s <- run_e_step_univariate(m, g, exp_m_offset, exp_g_offset, curr_m_perturbation, curr_g_perturbation, curr_pi)
    # M step
    m_step <- run_m_step_univariate(m, g, exp_m_offset, exp_g_offset, curr_Ti1s, n)
    # update log likelihood and estimates
    curr_log_lik <- m_step$log_lik
    curr_m_perturbation <- m_step$fit_m_perturbation; curr_g_perturbation <- m_step$fit_g_perturbation; curr_pi <- m_step$fit_pi
    # check for convergence
    curr_tol <- abs(curr_log_lik - prev_log_lik)/min(abs(curr_log_lik), abs(prev_log_lik))
    if (curr_tol < ep_tol) {
      converged <- TRUE
    } else {
      prev_log_lik <- curr_log_lik
      iteration <- iteration + 1L
      if (iteration >= max_it) break()
    }
  }
  out <- list(m_perturbation = curr_m_perturbation, g_perturbation = curr_g_perturbation, pi = curr_pi, log_lik = curr_log_lik, converged = converged, n_iterations = iteration)
  return(out)
}


run_e_step_univariate <- function(m, g, exp_m_offset, exp_g_offset, curr_m_perturbation, curr_g_perturbation, curr_pi) {
  quotient <- (log(1 - curr_pi) - log(curr_pi)) +
    (exp(curr_m_perturbation) - 1) * exp_m_offset - curr_m_perturbation * m +
    (exp(curr_g_perturbation) - 1) * exp_g_offset - curr_g_perturbation * g
  Ti1s <- 1/(exp(quotient) + 1)
  return(Ti1s)
}


run_m_step_univariate <- function(m, g, exp_m_offset, exp_g_offset, curr_Ti1s, n) {
  # estimate pi and compute pi log-likelihood
  s_curr_Ti1s <- sum(curr_Ti1s)
  fit_pi <- s_curr_Ti1s/n
  if (fit_pi >= 0.5) { # subtract by 1 to ensure label consistency
    s_curr_Ti1s <- n - s_curr_Ti1s
    fit_pi <- s_curr_Ti1s/n
  }
  pi_log_lik <- log(1 - fit_pi) * (n - s_curr_Ti1s) + log(fit_pi) * s_curr_Ti1s

  # fit the models for m and g
  fit_m <- fit_model_univariate(m, curr_Ti1s, exp_m_offset)
  fit_g <- fit_model_univariate(g, curr_Ti1s, exp_g_offset)
  # compute the overall model log-likelihood
  log_lik <- fit_m$log_lik + fit_g$log_lik + pi_log_lik
  # return fitted parameters, as well as log-likelihood
  out <- list(log_lik = log_lik, fit_pi = fit_pi,
              fit_m_perturbation = fit_m$fit_param, fit_g_perturbation = fit_g$fit_param)
  return(out)
}


fit_model_univariate <- function(v, curr_Ti1s, exp_offset) {
  p1 <- sum(curr_Ti1s * v)
  p2 <- sum(curr_Ti1s * exp_offset)
  fit_param <- log(p1/p2)
  log_lik <- fit_param * p1 - exp(fit_param) * p2
  return(list(fit_param = fit_param, log_lik = log_lik))
}
