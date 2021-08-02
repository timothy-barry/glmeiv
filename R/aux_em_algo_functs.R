# Consider depricating entire script.

run_em_algo_multiple_inits <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, return_best) {
  em_runs <- apply(X = initial_Ti1_matrix, MARGIN = 2, FUN = function(initial_Ti1s) {
    run_full_glmeiv_given_weights(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset)
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


#' Initialize weights using marginal mixtures
#'
#' Initializes weights for GLM-EIV by fitting marginal mixtures, then taking weighted average.
#'
#' @param m mRNA counts
#' @param g gRNA counts
#' @param m_fam family describing mRNA counts
#' @param g_fam family describing gRNA counts
#' @param m_offset optional offsets for m
#' @param g_offset optional offsets for g
#' @param lambda optional weight in weighted average; if not supplied, chosen adaptively
#'
#' @return initial weights for algorithm
initialize_weights_using_marginal_mixtures <- function(m, g, m_fam, g_fam, m_offset, g_offset, lambda = NULL) {
  m_weights <- get_marginal_mixture_weights(m, m_fam, m_offset)
  g_weights <- get_marginal_mixture_weights(g, g_fam, g_offset)
  if (is.null(lambda)) {
    d_m <- compute_mean_distance_from_half(m_weights)
    d_g <- compute_mean_distance_from_half(g_weights)
    d_sum <- d_m + d_g
    lambda <- if (d_sum < 1e-10) 1 else d_g/(d_sum)
  }
  out <- lambda * g_weights + (1 - lambda) * m_weights
  return(out)
}


#' Get marginal mixture weights
#'
#' Fit a marginal mixture model; get the (correctly oriented) weights
#'
#' @param v a vector (of mRNA or gRNA counts)
#' @param fam family object describing the counts
#' @param offset optional vector of linear offsets
#'
#' @return the EM weights
get_marginal_mixture_weights <- function(v, fam, offset) {
  flex_fit <- flexmix::flexmix(v ~ 1, k = 2,
                               model = flexmix::FLXglm(family = fam$flexmix_fam,
                                                       offset = offset))
  if (flex_fit@k == 1) {
    w <- rep(0.5, length(v))
  } else {
    w_matrix <- flex_fit@posterior$scaled
    w <- if (sum(w_matrix[,1]) <= sum(w_matrix[,2])) w_matrix[,1] else w_matrix[,2]
  }
  return(w)
}


#' Append noise to weights
#'
#' Adds Gaussian noise to an initial weight vector.
#'
#' @param w initial weight vector
#' @param n_rep number of noisy columns to append to w
#' @param sd standard deviation of noise
#'
#' @return initial weight matrix
append_noise_to_weights <- function(w, n_rep, sd) {
  initial_Ti1_matrix <- replicate(n = n_rep, {
    noise <- stats::rnorm(n = length(w), mean = 0, sd = sd)
    out <- w + noise
    out[out > 1] <- 1; out[out < 0] <- 0
    out
  }) %>% cbind(w, .)
  return(initial_Ti1_matrix)
}
