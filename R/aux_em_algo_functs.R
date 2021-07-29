run_em_algo_multiple_inits <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, return_best, ep_tol = 0.5 * 1e-4, max_it = 75) {
  em_runs <- apply(X = initial_Ti1_matrix, MARGIN = 2, FUN = function(initial_Ti1s) {
    run_em_algo_given_weights(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, ep_tol, max_it)
  })
  names(em_runs) <- NULL
  if (return_best) {
    out <- select_best_em_run(em_runs)
  } else {
    out <- em_runs
  }
  return(out)
}


#' Run EM algorithm -- mixture model initialization
#'
#' @param dat the count data; should have columns m and g.
#' @param g_fam family object used to model gRNA counts.
#' @param m_fam family object used to model mRNA counts.
#' @param covariate_matrix optional matrix of covariates
#' @param m_offset optional offsets for mRNA model
#' @param g_offset optional offsets for gRNA model
#' @param alpha returns a (1-alpha)% confidence interval
#' @param n_em_rep number of replicates of em algorithm
#' @param sd sd of noise to add to initial weights
#' @param save_membership_probs_mult save the posterior membership probabilities at this multiple
#' @param lambda mixing parameter for weighted average of weights; default NULL chooses mixing parameter adaptively.
#'
#' @return a tibble with columns parameter, target (fields estimate, std_error, p-value, confint lower, and confint higher), and value.
#' @export
#'
#' @examples
#' \dontrun{
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -1; g_perturbation <- 2
#' pi <- 0.2; n <- 2000; B <- 5; alpha <- 0.95; n_em_rep <- 3;
#' sd <- 0.15; lambda <- NULL
#' m_offset <- g_offset <- NULL
#' m_covariate_coefs <- g_covariate_coefs <- covariate_matrix <- NULL
#' dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
#' g_intercept, g_perturbation, pi, n, B, covariate_matrix,
#' m_covariate_coefs, g_covariate_coefs, m_offset, g_offset)
#' dat <- dat_list[[1]]
#' run_em_algo_mixture_init(dat, g_fam, m_fam, covariate_matrix,
#' m_offset, g_offset, alpha, n_em_rep, sd)
#' }
run_em_algo_mixture_init <- function(dat, g_fam, m_fam, covariate_matrix, m_offset, g_offset, alpha, n_em_rep, sd = 0.15, save_membership_probs_mult = 250, lambda = NULL) {
  m <- dat$m
  g <- dat$g
  # obtain initial weights
  w <- initialize_weights_using_marginal_mixtures(m = m, g = g, m_fam = m_fam, g_fam = g_fam, m_offset = m_offset, g_offset = g_offset, lambda = lambda)
  # obtain initial weight matrix by adding noise
  initial_Ti1_matrix <- append_noise_to_weights(w, n_em_rep - 1, sd)
  em_fit <- run_em_algo_multiple_inits(m = m, g = g, m_fam = m_fam, g_fam = g_fam, covariate_matrix = covariate_matrix,
                                       initial_Ti1_matrix = initial_Ti1_matrix, m_offset = m_offset,
                                       g_offset = g_offset, return_best = TRUE)
  # compute fit confidence metrics
  membership_prob_spread <- compute_mean_distance_from_half(em_fit$posterior_perturbation_probs)
  n_approx_1 <- sum(em_fit$posterior_perturbation_probs > 0.85)
  n_approx_0 <- sum(em_fit$posterior_perturbation_probs < 0.15)
  # do inference
  s <- run_inference_on_em_fit(em_fit, alpha) %>% dplyr::rename("parameter" = "variable")
  # output result
  meta_df <- tibble::tibble(parameter = "meta",
                            target = c("converged", "membership_probability_spread", "n_approx_0", "n_approx_1"),
                            value = c(em_fit$converged, membership_prob_spread, n_approx_0, n_approx_1))
  out <- rbind(tidyr::pivot_longer(s, cols = -parameter, names_to = "target"), meta_df)
  # if i is a multiple of 250, save the posterior membership probabilities
  i <- attr(dat, "i")
  if ((i - 1 + save_membership_probs_mult) %% save_membership_probs_mult == 0) {
    out <- rbind(out, data.frame(parameter = "meta",
                                 target = "membership_prob",
                                 value = em_fit$posterior_perturbation_probs))
  }
  return(out)
}


#' Run EM algo with fast init
#'
#' Runs the EM algorithm using the fast initialization strategy.
#' The assumption is that n is big and pi is small. It is also assumed
#' (for now) that the covariate matrix is non-null.
#'
#' NOTE: Maybe instead take precomputations, i.e. distillation offsets, as arguments.
#'
#'
#' @param dat a data frame containing the columns "m" and "g"
#' @param g_fam family object describing mRNA counts
#' @param m_fam family object describing gRNA counts
#' @param covariate_matrix a data frame storing the covariates
#' @param m_offset the vector of offsets for mRNA model
#' @param g_offset the vector of offsets for gRNA model
#' @param alpha confidence level of CIs
#' @param n_em_rep number of times to repeat EM algorithm on "reduced" data
#'
#' @return fitted model
#' @export
#' @examples
#' m_fam <- g_fam <- augment_family_object(poisson())
#' n <- 200000
#' lib_size <- rpois(n = n, lambda = 10000)
#' m_offset <- g_offset <- log(lib_size)
#' pi <- 0.005
#' m_intercept <- log(0.05)
#' m_perturbation <- log(0.75)
#' g_intercept <- log(0.025)
#' g_perturbation <- log(1.25)
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
#' fit <- run_em_algo_fast_init(m, g, m_fam, g_fam, covariate_matrix, m_offset, g_offset)
run_em_algo_fast_init <- function(m, g, m_fam, g_fam, covariate_matrix, m_offset, g_offset, n_em_rep = 15, pi_guess_range = c(1e-5, 0.02), m_perturbation_guess_range = log(c(0.1, 1.5)), g_perturbation_guess_range = log(c(0.5, 10)), alpha = 0.95) {
  # run the mRNA and gRNA precomputations
  fit_m_precomp <- stats::glm(formula = m ~ ., family = m_fam, data = covariate_matrix, offset = m_offset)
  fit_m_precomp_coefs <- coef(fit_m_precomp)
  fit_g_precomp <- stats::glm(formula = g ~ ., family = g_fam, data = covariate_matrix, offset = g_offset)
  fit_g_precomp_coefs <- coef(fit_g_precomp)

  # obtain the fitted values
  fitted_vals_m_precomp <- as.numeric(stats::fitted.values(fit_m_precomp))
  fitted_vals_g_precomp <- as.numeric(stats::fitted.values(fit_g_precomp))

  # obtain random starting points for reduced GLM-EIV model
  set.seed(4)
  pi_guess <- runif(n = n_em_rep, min = pi_guess_range[1], max = pi_guess_range[2])
  m_perturbation_guess <- runif(n = n_em_rep, min = m_perturbation_guess_range[1], max = m_perturbation_guess_range[2])
  g_perturbation_guess <- runif(n = n_em_rep, min = g_perturbation_guess_range[1], max = g_perturbation_guess_range[2])

  # fit the reduced GLM-EIV model n_em_rep times
  reduced_fits <- lapply(seq(1L, n_em_rep), function(i) {
    run_univariate_poisson_em_algo(m = m, g = g, exp_m_offset = fitted_vals_m_precomp, exp_g_offset = fitted_vals_g_precomp,
                                   m_pert_guess = m_perturbation_guess[i], g_pert_guess = g_perturbation_guess[i], pi_guess = pi_guess[i])
  })

  # among fits that converged, select the one with greatest log-likelihood
  converged <- sapply(reduced_fits, function(fit) fit$converged)
  reduced_fit <- select_best_em_run(reduced_fits[converged])

  # obtain membership probabilities to initialize EM algo
  technical_factors <- colnames(covariate_matrix)
  initial_Ti1s <- run_e_step_pilot(m = m,
                                   g = g,
                                   m_fam = m_fam,
                                   g_fam = g_fam,
                                   pi_guess = reduced_fit$pi,
                                   m_intercept_guess = fit_m_precomp_coefs[["(Intercept)"]],
                                   m_perturbation_guess = reduced_fit$m_perturbation,
                                   m_covariate_coefs_guess = fit_m_precomp_coefs[[technical_factors]],
                                   g_intercept_guess = fit_g_precomp_coefs[["(Intercept)"]],
                                   g_perturbation_guess = reduced_fit$g_perturbation,
                                   g_covariate_coefs_guess = fit_g_precomp_coefs[[technical_factors]],
                                   covariate_matrix = covariate_matrix,
                                   m_offset = m_offset,
                                   g_offset = g_offset)

  # run em algo with initial weights
  fit_em <- run_em_algo_given_weights(m = m, g = g, m_fam = m_fam, g_fam = g_fam, covariate_matrix = covariate_matrix,
                                      initial_Ti1s = initial_Ti1s$Ti1s, m_offset = m_offset, g_offset = g_offset, prev_log_lik = initial_Ti1s$log_lik)
  out <- list(run_inference_on_em_fit(fit_em, alpha = alpha), fit_em)
  return(out)
}
