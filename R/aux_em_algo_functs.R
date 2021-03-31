#' @rdname run_em_algo_given_init
#' @export
run_em_algo_multiple_inits <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, ep_tol = 0.1, max_it = 50) {
  em_runs <- apply(X = initial_Ti1_matrix, MARGIN = 2, FUN = function(initial_Ti1s) {
    run_em_algo_given_init(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, ep_tol = ep_tol, max_it = max_it)
  })
}


#' Select best EM run
#'
#' Returns the best run, defined as the run with positive perturbation coefficient for G and
#' greatest log-likelihood.
#'
#' @param em_runs a list of outputs of run_em_algo_given_init
#'
#' @return the best run of the list
#' @export
#' @examples
#' data <- get_quick_simulated_data()
select_best_em_run <- function(em_runs) {
  runwise_stats <- sapply(em_runs, function(run)
    c(log_lik = run$log_lik, pert_g = stats::coef(run$fit_g)[["perturbation"]])) %>%
    rbind(run_id = seq(1, length(em_runs)))
  pos_g_pert <- runwise_stats[,runwise_stats["pert_g",] >= 0]
  best_idx <- pos_g_pert[["run_id", which.max(pos_g_pert["log_lik",])]] %>% as.integer()
  return(em_runs[[best_idx]])
}
