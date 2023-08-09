#' Run unimodal mixture method in simulatr
#' @export
#' @examples
#' m_fam <- g_fam <- augment_family_object(poisson())
#' n <- 50000
#' lib_size <- rpois(n = n, lambda = 5000)
#' m_offset <- g_offset <- log(lib_size)
#' pi <- 0.01
#' m_intercept <- log(0.01)
#' m_perturbation <- log(0.5)
#' g_intercept <- log(0.005)
#' g_perturbation <- log(2.5)
#' covariate_matrix <- data.frame(batch = rbinom(n = n, size = 1, prob = 0.5))
#' m_covariate_coefs <- log(0.9)
#' g_covariate_coefs <- log(1.1)
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept,
#' m_perturbation = m_perturbation, g_fam = g_fam, g_intercept = g_intercept,
#' g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = covariate_matrix, m_covariate_coefs = m_covariate_coefs,
#' g_covariate_coefs = g_covariate_coefs, m_offset = m_offset, g_offset = g_offset)[[1]]
#' m <- dat$m; g <- dat$g; p <- dat$p
#' fit <- run_umimodal_mixture_method_simulatr(dat, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset)
run_umimodal_mixture_method_simulatr <- function(dat, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha = 0.95, n_em_rep = 15, pi_guess_range = c(1e-5, 0.03), g_perturbation_guess_range = log(c(0.5, 10))) {
  set.seed(4)
  if (!is(m_fam, "family")) m_fam <- m_fam[[1]]
  if (!is(g_fam, "family")) g_fam <- g_fam[[1]]
  # pull g_fam and m_fam from dat, if available
  if ("m_precomp" %in% names(attributes(dat))) {
    m_precomp <- attr(dat, "m_precomp"); m_fam <- m_precomp$fam
  }
  if ("g_precomp" %in% names(attributes(dat))) {
    g_precomp <- attr(dat, "g_precomp"); g_fam <- g_precomp$fam
  } else {
    g_precomp <- run_glmeiv_precomputation(y = g, covariate_matrix = covariate_matrix, offset = g_offset, fam = g_fam)
    g_fam <- g_precomp$fam
  }

  time <- system.time({
    # step 0: get random starting guesses and fitted grna regression values
    guesses <- lapply(list(pi = pi_guess_range,
                           g_perturbation = g_perturbation_guess_range), function(r) {
                             stats::runif(n = n_em_rep, min = r[1], max = r[2])})
    g_fitted <- g_fam$linkfun(g_precomp$fitted_values)

    # step 1: run the reduced glmeiv using the gRNA modality only
    reduced_fits <- lapply(seq(1L, n_em_rep), function(i) {
      run_reduced_em_algo(m = NULL, g = g, m_fitted = NULL, g_fitted = g_fitted,
                          m_pert_guess = NULL,
                          g_pert_guess = guesses$g_perturbation[i],
                          pi_guess = guesses$pi[i], m_fam = NULL, g_fam = g_fam, use_mrna_modality = FALSE)
    })

    # step 2: obtain the best run according to log-likelihood
    log_liks <- sapply(X = reduced_fits, FUN = function(l) l[["log_lik"]])
    best_run <- which.max(log_liks)
    best_fit <- reduced_fits[[best_run]]

    # step 3: run the full unimodal grna mixture model
    fit <- run_full_glmeiv_given_pilot_params(m = NULL, g = g, m_fam = NULL, g_fam = g_fam,
                                              pi_guess = best_fit$pi, m_intercept_guess = NULL,
                                              m_perturbation_guess = NULL, m_covariate_coefs_guess = NULL,
                                              g_intercept_guess = g_precomp$fitted_intercept, g_perturbation_guess = best_fit$g_perturbation,
                                              g_covariate_coefs_guess = g_precomp$covariate_coefs, covariate_matrix = covariate_matrix,
                                              m_offset = m_offset, g_offset = g_offset, use_mrna_modality = FALSE)

    # step 4: extract the clusters (phat) and fit the GLM
    phat <- as.integer(fit$posterior_perturbation_probs >= 0.5)

    # step 5: call the thresholding method simulatr helper
    m <- dat$m
    out <- thresholding_method_simulatr_helper(phat, m, pi, covariate_matrix, m_fam, m_offset, alpha)
  })[["elapsed"]]
  out <- out %>% dplyr::add_row(.,parameter = "meta", target = "time", value = time)
  return(out)
}
