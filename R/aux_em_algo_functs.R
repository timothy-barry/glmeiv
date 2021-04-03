# This script contains "auxiliary" EM algorithm functions. The functions here are more experimental than those in other R files in this package.

#' @rdname run_em_algo_given_init
#' @export
run_em_algo_multiple_inits <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, ep_tol = 0.1, max_it = 50) {
  em_runs <- apply(X = initial_Ti1_matrix, MARGIN = 2, FUN = function(initial_Ti1s) {
    run_em_algo_given_init(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, ep_tol = ep_tol, max_it = max_it)
  })
  names(em_runs) <- NULL
  return(em_runs)
}


#' Select best EM run
#'
#' Returns the best run, defined as the run with positive perturbation coefficient for G and
#' greatest log-likelihood.
#'
#' @param em_runs a list of outputs of run_em_algo_given_init
#'
#' @return the best run of the list
select_best_em_run <- function(em_runs) {
  # select the run with greatest likelihood
  log_liks <- sapply(em_runs, function(run) run$log_lik)
  return(em_runs[[which.max(log_liks)]])
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
#'
#' @return the best EM run
#' @export
#'
#' @examples
#' dat <- get_quick_simulated_data(5000)
#' em_coefs <- run_glmeiv_known_p(m = dat$m, g = dat$g, m_fam = dat$m_fam, g_fam = dat$m_fam,
#' covariate_matrix = dat$covariate_matrix, p = dat$p, m_offset = NULL, g_offset = NULL)
#' # dat$m_coef; dat$g_coef; dat$pi
run_glmeiv_known_p <- function(m, g, m_fam, g_fam, covariate_matrix, p, m_offset = NULL, g_offset = NULL, n_runs = 4, p_flip = 0.15, ep_tol = 0.1, max_it = 50) {
  initial_Ti1_matrix <- cbind(p, replicate(n_runs - 1, expr = flip_weights(p, p_flip)))
  em_runs <- run_em_algo_multiple_inits(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, ep_tol = ep_tol, max_it = max_it)
  best_run <- select_best_em_run(em_runs)
  out <- run_inference_on_em_fit(best_run)
  return(out)
}
