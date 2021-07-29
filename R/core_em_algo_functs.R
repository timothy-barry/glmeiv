#' Run EM algo given initialization weights
#'
#' Runs GLM-EIV algo with M-step first. Use `run_em_algo_given_pilot` to GLM-EIV algo with E-step first.
#'
#' @param m observed vector m
#' @param g observed vector g
#' @param m_fam augmented family of m
#' @param g_fam augmented family of g
#' @param covariate_matrix the matrix of covariates; set to NULL if there are no covariates.
#' @param initial_Ti1s the initial vector of membership probabilities
#' @param m_offset offsets for GLM for M
#' @param g_offset offsets for GLM for G
#' @param ep_tol (optional) EM convergence threshold
#' @param max_it  (optional) maximum number of EM iterations
#'
#' @return a list containing the following: fit object of M GLM, fit object of G GLM,
#' fit for pi, number of iterations,full log-likelihood of final model
#' @export
#' @name run_em_algo_given_init
#' @examples
#' m_fam <- g_fam <- augment_family_object(poisson())
#' n <- 5000
#' lib_size <- rpois(n = n, lambda = 5000)
#' m_offset <- g_offset <- log(lib_size)
#' pi <- 0.1
#' m_intercept <- log(0.05)
#' m_perturbation <- log(0.8)
#' g_intercept <- log(0.025)
#' g_perturbation <- log(1.2)
#' covariate_matrix <- data.frame(batch = rbinom(n = n, size = 1, prob = 0.5))
#' m_covariate_coefs <- log(0.9)
#' g_covariate_coefs <- log(1.1)
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept,
#' m_perturbation = m_perturbation, g_fam = g_fam, g_intercept = g_intercept,
#' g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = covariate_matrix, m_covariate_coefs = m_covariate_coefs,
#' g_covariate_coefs = g_covariate_coefs, m_offset = m_offset, g_offset = g_offset)[[1]]
#' pi_guess <- 0.05
#' m_intercept_guess <- log(0.07)
#' m_perturbation_guess <- log(0.7)
#' g_intercept_guess <- log(0.02)
#' g_perturbation_guess <- log(1.4)
#' m_covariate_coefs_guess <- log(0.8)
#' g_covariate_coefs_guess <- log(1.2)
#' m <- dat$m
#' g <- dat$g
#' # obtain initial membership probabilities (i.e., run E step) using pilot estimates
#' initial_Ti1s <- run_e_step_pilot(m, g, m_fam, g_fam, pi_guess,
#' m_intercept_guess, m_perturbation_guess, m_covariate_coefs_guess,
#' g_intercept_guess, g_perturbation_guess, g_covariate_coefs_guess,
#' covariate_matrix, m_offset, g_offset)
#' # run em algo
#' fit <- run_em_algo_given_weights(m, g, m_fam, g_fam, covariate_matrix,
#' initial_Ti1s, m_offset, g_offset)
run_em_algo_given_weights <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, prev_log_lik = -Inf, ep_tol = 1e-4, max_it = 75) {
  # augment family objects, if necessary
  if (is.null(m_fam$augmented)) m_fam <- augment_family_object(m_fam)
  if (is.null(g_fam$augmented)) g_fam <- augment_family_object(g_fam)

  # verify column names ok
  check_col_names(covariate_matrix)

  # define some basic quantities
  n <- length(m)
  iteration <- 1L
  converged <- FALSE
  curr_Ti1s <- initial_Ti1s
  prev_log_lik <- prev_log_lik
  log_liks <- numeric()

  # define augmented responses, offsets, and covariate matrix
  augmented_inputs <- augment_inputs(covariate_matrix, m, g, m_offset, g_offset, n)

  # iterate through E and M steps until convergence
  while (!converged) {
    # m step
    m_step <- run_m_step(curr_Ti1s,
                         augmented_inputs$m_augmented, m_fam, augmented_inputs$m_offset_augmented,
                         augmented_inputs$g_augmented, g_fam, augmented_inputs$g_offset_augmented,
                         augmented_inputs$Xtilde_augmented, n)
    curr_log_lik <- m_step$curr_log_lik
    log_liks <- c(log_liks, curr_log_lik)
    curr_tol <- abs(curr_log_lik - prev_log_lik)/min(abs(curr_log_lik), abs(prev_log_lik))
    if (curr_tol < ep_tol) {
      # convergence acheived
      converged <- TRUE
    } else {
      # e step
      curr_Ti1s <- run_e_step(m_step, m, m_fam, g, g_fam, n)
      prev_log_lik <- curr_log_lik
      iteration <- iteration + 1L
      # check iteration limit
      if (iteration >= max_it) {
        break()
      }
    }
  }
  out <- list(fit_m = m_step$fit_m, fit_g = m_step$fit_g, fit_pi = m_step$fit_pi, n_iterations = iteration, log_liks = log_liks, log_lik = curr_log_lik, converged = converged, n = n, posterior_perturbation_probs = m_step$posterior_perturbation_probs)
  return(out)
}


##################
# helper functions
##################
augment_inputs <- function(covariate_matrix, m, g, m_offset, g_offset, n) {
  if (is.null(covariate_matrix)) {
    Xtilde_augmented <- data.frame(perturbation = c(rep(0, n), rep(1, n)))
  } else {
    Xtilde_0 <- dplyr::mutate(covariate_matrix, perturbation = 0)
    Xtilde_1 <- dplyr::mutate(covariate_matrix, perturbation = 1)
    Xtilde_augmented <- rbind(Xtilde_0, Xtilde_1) %>% dplyr::select(perturbation, everything())
  }
  m_augmented <- c(m, m)
  g_augmented <- c(g, g)
  m_offset_augmented <- if (!is.null(m_offset)) c(m_offset, m_offset) else NULL
  g_offset_augmented <- if (!is.null(g_offset)) c(g_offset, g_offset) else NULL
  out <- list(Xtilde_augmented = Xtilde_augmented,
              m_augmented = m_augmented,
              g_augmented = g_augmented,
              m_offset_augmented = m_offset_augmented,
              g_offset_augmented = g_offset_augmented)
  return(out)
}


run_m_step <- function(curr_Ti1s, m_augmented, m_fam, m_offset_augmented, g_augmented, g_fam, g_offset_augmented, Xtilde_augmented, n) {
  # curr_Ti1s[curr_Ti1s < 1e-4] <- 0 # induce weight sparsity
  weights <- c(1 - curr_Ti1s, curr_Ti1s)
  s_curr_Ti1s <- sum(curr_Ti1s)
  fit_pi <- s_curr_Ti1s/n
  if (fit_pi >= 0.5) { # subtract by 1 to ensure label consistency
    s_curr_Ti1s <- n - s_curr_Ti1s
    fit_pi <- s_curr_Ti1s/n
  }
  pi_log_lik <- log(1 - fit_pi) * (n - s_curr_Ti1s) + log(fit_pi) * s_curr_Ti1s

  # fit models for m and g
  m_form <- stats::formula("m_augmented ~ .")
  fit_m <- stats::glm(formula = m_form, data = Xtilde_augmented, family = m_fam,
                      weights = weights, offset = m_offset_augmented)

  g_form <-  stats::formula("g_augmented ~ .")
  fit_g <- stats::glm(formula = g_form, data = Xtilde_augmented, family = g_fam,
                      weights = weights, offset = g_offset_augmented)

  # compute the log-likelihoods
  m_log_lik <- m_fam$get_log_lik(fit_m)
  g_log_lik <- g_fam$get_log_lik(fit_g)

  curr_log_lik <- m_log_lik + g_log_lik + pi_log_lik

  # return list of fitted models, as well as current log-likelihood
  out <- list(fit_m = fit_m, fit_g = fit_g, fit_pi = fit_pi, curr_log_lik = curr_log_lik, posterior_perturbation_probs = curr_Ti1s)
  return(out)
}


update_membership_probs <- function(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi) {
  quotient <- log(1 - fit_pi) - log(fit_pi) + m_fam$d_log_py(m, m_mus_pert0, m_mus_pert1) + g_fam$d_log_py(g, g_mus_pert0, g_mus_pert1)
  out <- 1/(exp(quotient) + 1)
  return(out)
}


run_e_step <- function(m_step, m, m_fam, g, g_fam, n) {
  compute_conditional_means <- function(fit, fam, n) {
    mus <- fam$linkinv(as.numeric(fit$linear.predictors))
    mus_pert0 <- mus[seq(1, n)]
    mus_pert1 <- mus[seq(n + 1, 2 * n)]
    out <- list(mus_pert0 = mus_pert0, mus_pert1 = mus_pert1)
  }
  # compute conditional means
  m_mus <- compute_conditional_means(m_step$fit_m, m_fam, n)
  g_mus <- compute_conditional_means(m_step$fit_g, g_fam, n)
  # define all relevant variables
  fit_pi <- m_step$fit_pi
  m_mus_pert0 <- m_mus$mus_pert0; m_mus_pert1 <- m_mus$mus_pert1
  g_mus_pert0 <- g_mus$mus_pert0; g_mus_pert1 <- g_mus$mus_pert1
  # compute membership probabilities
  Ti1s <- update_membership_probs(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi)
  return(Ti1s)
}


run_e_step_pilot <- function(m, g, m_fam, g_fam, pi_guess, m_intercept_guess, m_perturbation_guess, m_covariate_coefs_guess, g_intercept_guess, g_perturbation_guess, g_covariate_coefs_guess, covariate_matrix, m_offset, g_offset) {
  # compute the conditional means
  m_conditional_means <- compute_theoretical_conditional_means(intercept = m_intercept_guess,
                                                               perturbation_coef = m_perturbation_guess,
                                                               fam = m_fam,
                                                               covariate_matrix = covariate_matrix,
                                                               covariate_coefs = m_covariate_coefs_guess,
                                                               offset = m_offset)
  g_conditional_means <- compute_theoretical_conditional_means(intercept = g_intercept_guess,
                                                               perturbation_coef = g_perturbation_guess,
                                                               fam = g_fam,
                                                               covariate_matrix = covariate_matrix,
                                                               covariate_coefs = g_covariate_coefs_guess,
                                                               offset = g_offset)
  # assign to variables for convenience
  m_mus_pert0 <- m_conditional_means$mu0; m_mus_pert1 <- m_conditional_means$mu1
  g_mus_pert0 <- g_conditional_means$mu0; g_mus_pert1 <- g_conditional_means$mu1
  # compute membership probabilities
  Ti1s <- update_membership_probs(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, pi_guess)
  # compute the log-likelihood
  log_lik <- compute_weighted_log_lik(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, pi_guess, Ti1s)
  return(list(Ti1s = Ti1s, log_lik = log_lik))
}


compute_weighted_log_lik <- function(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi, Ti1s) {
  # weighted pi log-likelihood; assumes that the Ti1s are correct, i.e., sum(Ti1s)/n < 0.5.
  s_curr_Ti1s <- sum(Ti1s)
  pi_log_lik <- log(1 - fit_pi) * (n - s_curr_Ti1s) + log(fit_pi) * s_curr_Ti1s
  # mRNA log likelihood
  mRNA_log_lik <- m_fam$weighted_log_lik(y = m, mu_0 = m_mus_pert0, mu_1 = m_mus_pert1, Ti1s = Ti1s)
  # gRNA log likelihood
  gRNA_log_lik <- g_fam$weighted_log_lik(y = g, mu_0 = g_mus_pert0, mu_1 = g_mus_pert1, Ti1s = Ti1s)
  # log likelihood
  log_lik <- pi_log_lik + mRNA_log_lik + gRNA_log_lik
  return(log_lik)
}
