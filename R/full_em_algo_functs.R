#' Run E step
#'
#' Runs the E step of the GLM-EIV EM algorithm. Also, returns the updated value of pi (thereby running a portion of the subsequent M step).
#'
#' Ensures that the value of pi in the next M iteration will be less than 0.5. This core function is called in all variants of GLM-EIV (weights/estimates first, full/reduced, etc.).
#'
#' @param m_fam family object describing mRNA counts
#' @param g_fam family object describing gRNA counts
#' @param m mRNA counts
#' @param g gRNA counts
#' @param m_mus_pert0 fitted (or hypothesized) means mu^m(0)
#' @param m_mus_pert1 fitted (or hypothesized) means mu^m(1)
#' @param g_mus_pert0 fitted (or hypothesized) means mu^g(0)
#' @param g_mus_pert1 fitted (or hypothesized) means mu^g(1)
#' @param fit_pi fitted (or hypothesized) value for pi
#'
#' @return a list containing (i) the membership probabilities (Ti1s), (ii) the model log-likelihood, and (iii) the new value of pi (computed ahead of the subsequent M step for convenience).
run_e_step <- function(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi) {
  # first, compute log-likelihood
  p0 <- exp(log(1 - fit_pi) + m_fam$log_py_given_mu(m, m_mus_pert0) + g_fam$log_py_given_mu(g, g_mus_pert0))
  p1 <- exp(log(fit_pi) + m_fam$log_py_given_mu(m, m_mus_pert1) + g_fam$log_py_given_mu(g, g_mus_pert1))
  log_lik <- sum(log(p0 + p1))

  # second, compute membership probabilities
  quotient <- log(1 - fit_pi) - log(fit_pi) + m_fam$d_log_py(m, m_mus_pert0, m_mus_pert1) + g_fam$d_log_py(g, g_mus_pert0, g_mus_pert1)
  Ti1s <- 1/(exp(quotient) + 1)
  # estimate new_pi and ensure new_pi is less than 0.5
  new_pi <- sum(Ti1s)/length(Ti1s)
  if (new_pi > 0.5) {
    Ti1s <- 1 - Ti1s
    new_pi <- 1 - new_pi
  }
  # return quantities
  return(list(Ti1s = Ti1s, log_lik = log_lik, new_pi = new_pi))
}


#' Run full GLM-EIV given weights
#'
#' Runs the full GLM-EIV model given starting weights. Starting
#'
#' @param m mRNA counts
#' @param g gRNA counts
#' @param m_fam family describing m
#' @param g_fam family describing g
#' @param covariate_matrix the matrix of covariates; NULL if there are no covariates.
#' @param initial_Ti1s starting membership probabilities; these probabilities should be such that sum(initial_Ti1s)/length(initial_Ti1s) < 0.5.
#' @param prev_log_lik optional starting log-likelihood value (useful if an E step has been called ahead of running this function).
#' @param m_offset offsets for m
#' @param g_offset offsets for g
#' @param ep_tol (optional) EM convergence threshold
#' @param max_it  (optional) maximum number of EM iterations
#'
#' @return
#' @export
#'
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
#' m <- dat$m
#' g <- dat$g
#' initial_Ti1s <- runif(n)
#' fit <- run_full_glmeiv_given_weights(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset)
run_full_glmeiv_given_weights <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, prev_log_lik = -Inf, ep_tol = 1e-5, max_it = 75) {
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
  log_liks <- if (prev_log_lik == -Inf) numeric() else prev_log_lik
  augmented_inputs <- augment_inputs(covariate_matrix, m, g, m_offset, g_offset, n)

  # iterate through M and E steps (in that order) until convergence
  while (!converged) {
    # M step
    m_step <- run_m_step_full(curr_Ti1s, augmented_inputs$m_augmented, m_fam,
                              augmented_inputs$m_offset_augmented,
                              augmented_inputs$g_augmented, g_fam,
                              augmented_inputs$g_offset_augmented,
                              augmented_inputs$Xtilde_augmented, n)
    # E step
    e_step <- run_e_step(m_fam, g_fam, m, g, m_step$m_mus_pert0, m_step$m_mus_pert1,
                         m_step$g_mus_pert0, m_step$g_mus_pert1, m_step$fit_pi)
    # append current log-likelihood, check for convergence
    curr_Ti1s <- e_step$Ti1s
    curr_log_lik <- e_step$log_lik
    log_liks <- c(log_liks, curr_log_lik)
    tol <- compute_tolerance(curr_log_lik, prev_log_lik)
    if (tol < ep_tol) {
      converged <- TRUE
    } else {
      prev_log_lik <- curr_log_lik
      iteration <- iteration + 1L
      if (iteration >= max_it) break()
    }
  }
  out <- list(fit_m = m_step$fit_m, fit_g = m_step$fit_g, fit_pi = m_step$fit_pi,
              n_iterations = iteration, log_liks = log_liks, log_lik = curr_log_lik,
              converged = converged, n = n, posterior_perturbation_probs = curr_Ti1s)
  return(out)
}


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


run_m_step_full <- function(curr_Ti1s, m_augmented, m_fam, m_offset_augmented, g_augmented, g_fam, g_offset_augmented, Xtilde_augmented, n) {
  fit_pi <- sum(curr_Ti1s)/n
  # fit models for m and g using weights
  weights <- c(1 - curr_Ti1s, curr_Ti1s)
  fit_m <- stats::glm(formula = stats::formula("m_augmented ~ ."), data = Xtilde_augmented,
                      family = m_fam, weights = weights, offset = m_offset_augmented)
  fit_g <- stats::glm(formula = stats::formula("g_augmented ~ ."), data = Xtilde_augmented,
                      family = g_fam, weights = weights, offset = g_offset_augmented)
  # extract fitted means
  fitted_means_m <- extract_glm_fitted_means(fit_m, m_fam, n)
  fitted_means_g <- extract_glm_fitted_means(fit_g, g_fam, n)
  # return list fitted means, as well as the fitted GLM objects themselves
  out <- list(fit_m = fit_m, fit_g = fit_g, fit_pi = fit_pi,
              m_mus_pert0 = fitted_means_m$mus_pert0, m_mus_pert1 = fitted_means_m$mus_pert1,
              g_mus_pert0 = fitted_means_g$mus_pert0, g_mus_pert1 = fitted_means_g$mus_pert1)
  return(out)
}


# extract the fitted means
extract_glm_fitted_means <- function(fit, fam, n) {
  mus <- fam$linkinv(as.numeric(fit$linear.predictors))
  mus_pert0 <- mus[seq(1, n)]
  mus_pert1 <- mus[seq(n + 1, 2 * n)]
  out <- list(mus_pert0 = mus_pert0, mus_pert1 = mus_pert1)
  return(out)
}


compute_tolerance <- function(curr_log_lik, prev_log_lik) {
  if (curr_log_lik == -Inf) {
    tol <- Inf
  } else {
    tol <- abs(curr_log_lik - prev_log_lik)/min(abs(curr_log_lik), abs(prev_log_lik))
  }
  return(tol)
}


#' Run full GLM-EIV given fitted means
#'
#' @return
#' @export
#' @inheritParams run_full_glmeiv_given_weights
#' @inheritParams run_e_step
run_full_glmeiv_given_fitted_means <- function(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi, covariate_matrix, m_offset, g_offset, ep_tol = 1e-4, max_it = 75) {
  e_step <- run_e_step(m_fam = m_fam, g_fam = g_fam, m = m, g = g,
                       m_mus_pert0 = m_mus_pert0, m_mus_pert1 = m_mus_pert1,
                       g_mus_pert0 = g_mus_pert0, g_mus_pert1 = g_mus_pert1, fit_pi = fit_pi)
  out <- run_full_glmeiv_given_weights(m = m, g = g, m_fam = m_fam, g_fam = g_fam,
                                covariate_matrix = covariate_matrix, initial_Ti1s = e_step$Ti1s,
                                m_offset = m_offset, g_offset = g_offset, prev_log_lik = e_step$log_lik,
                                ep_tol = ep_tol, max_it = max_it)
  return(out)
}


#' Run full GLM-EIV given pilot estiamtes
#'
#' @param pi_guess pilot guess for pi
#' @param m_intercept_guess pilot guess for m_intercept
#' @param m_perturbation_guess pilot guess for m_perturbation
#' @param m_covariate_coefs_guess pilot guess for m_covariate_coefficients
#' @param g_intercept_guess pilot guess for g_intercept
#' @param g_perturbation_guess pilot guess for g_perturbation
#' @param g_covariate_coefs_guess pilot guess for g_covariate_coefs
#' @return
#' @export
#'
#' @inheritParams run_full_glmeiv_given_weights
#' @examples
#' set.seed(4)
#' m_fam <- g_fam <- augment_family_object(poisson())
#' n <- 5000
#' lib_size <- rpois(n = n, lambda = 5000)
#' m_offsets <- g_offsets <- log(lib_size)
#' pi <- 0.3
#' m_intercept <- log(0.05)
#' m_perturbation <- log(0.75)
#' g_intercept <- log(0.025)
#' g_perturbation <- log(1.4)
#' covariate_matrix <- data.frame(batch = rbinom(n = n, size = 1, prob = 0.5))
#' m_covariate_coefs <- log(0.9)
#' g_covariate_coefs <- log(1.1)
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept,
#' m_perturbation = m_perturbation, g_fam = g_fam, g_intercept = g_intercept,
#' g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = covariate_matrix, m_covariate_coefs = m_covariate_coefs,
#' g_covariate_coefs = g_covariate_coefs, m_offset = m_offsets, g_offset = g_offsets)[[1]]
#' m <- dat$m
#' g <- dat$g
#' fit <- run_full_glmeiv_given_pilot_params(m = m, g = g, m_fam = m_fam, g_fam = g_fam,
#' pi_guess = 0.15, m_intercept_guess = log(0.1), m_perturbation_guess = log(1),
#' m_covariate_coefs_guess = log(1.4), g_intercept_guess = log(0.05),
#' g_perturbation_guess = log(1.4), g_covariate_coefs_guess = log(1.2),
#' covariate_matrix = covariate_matrix, m_offset = m_offsets, g_offset = g_offsets)
run_full_glmeiv_given_pilot_params <- function(m, g, m_fam, g_fam, pi_guess, m_intercept_guess, m_perturbation_guess, m_covariate_coefs_guess, g_intercept_guess, g_perturbation_guess, g_covariate_coefs_guess, covariate_matrix, m_offset, g_offset, ep_tol = 1e-4, max_it = 75) {
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
  # run glmeiv given fitted means
  out <- run_full_glmeiv_given_fitted_means(m_fam = m_fam, g_fam = g_fam, m = m, g = g,
                                            m_mus_pert0 = m_mus_pert0, m_mus_pert1 = m_mus_pert1,
                                            g_mus_pert0 = g_mus_pert0, g_mus_pert1 = g_mus_pert1,
                                            fit_pi = pi_guess, covariate_matrix = covariate_matrix,
                                            m_offset = m_offset, g_offset = g_offset, ep_tol = ep_tol, max_it = max_it)
  return(out)
}


compute_mean_distance_from_half <- function(v) {
  n <- length(v)
  sum(abs(v - 0.5))/n
}
