#' @rdname run_em_algo_given_init
#' @export
run_em_algo_multiple_inits <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, return_best, ep_tol = 0.5 * 1e-4, max_it = 75) {
  em_runs <- apply(X = initial_Ti1_matrix, MARGIN = 2, FUN = function(initial_Ti1s) {
    run_em_algo_given_init(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, ep_tol, max_it)
  })
  names(em_runs) <- NULL
  if (return_best) {
    out <- select_best_em_run(em_runs)
  } else {
    out <- em_runs
  }
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
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 1
#' pi <- 0.2; n <- 1000; B <- 500; alpha <- 0.95; n_em_rep <- 5; p_flip <- 0.01
#' m_offset <- g_offset <- NULL
#' m_covariate_coefs <- g_covariate_coefs <- covariate_matrix <- NULL
#' dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
#' g_intercept, g_perturbation, pi, n, B, covariate_matrix,
#' m_covariate_coefs, g_covariate_coefs, m_offset, g_offset)
#' run_em_algo_simulatr_optimal_thresh(dat_list, g_intercept, g_perturbation,
#' g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset,
#' alpha, n_em_rep, p_flip)
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
      dplyr::mutate(run_id = i)
  })
  return(do.call(rbind, res_list))
}
