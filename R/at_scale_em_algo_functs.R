#' Run GLM-EIV given precomputations
#'
#' Runs a full GLM-EIV given mRNA and gRNA precomputations, as outputted by `run_glmeiv_precomputation` function.
#'
#' @param m_precomp mRNA precomputation list, as outputted by `run_glmeiv_precomputation`
#' @param g_precomp gRNA precomputation list, as outputted by `run_glmeiv_precomputation`
#' @param fam_str string describing family; currently "poisson" or "Negative Binomial" supported
#' @param n_em_rep number of EM reps to run on reduced GLM-EIV model
#' @param pi_guess_range range of values over which to sample pi
#' @param m_perturbation_guess_range range of values over which to sample m_perturbation
#' @param g_perturbation_guess_range range of values over which to sample g_perturbation
#' @inheritParams run_full_glmeiv_given_weights
#'
#' @return a fitted GLM-EIV model
#' @export
run_glmeiv_given_precomputations <- function(m, g, m_precomp, g_precomp, fam_str, covariate_matrix, m_offset, g_offset, n_em_rep, pi_guess_range, m_perturbation_guess_range, g_perturbation_guess_range) {
  set.seed(4)
  # set family objects
  if (fam_str == "poisson") {
    m_fam <- g_fam <- augment_family_object(poisson())
  } else if (fam_str == "Negative Binomial") {
    m_fam <- augment_family_object(MASS::negative.binomial(m_precomp$theta))
    g_fam <- augment_family_object(MASS::negative.binomial(g_precomp$theta))
  }

  # get random starting guesses
  guesses <- lapply(list(pi = pi_guess_range,
                         m_perturbation = m_perturbation_guess_range,
                         g_perturbation = g_perturbation_guess_range), function(r) {
                           stats::runif(n = n_em_rep, min = r[1], max = r[2])
                         })

  # fit the reduced GLM-EIV model over the random starting vectors
  exp_m_offset <- m_precomp$fitted_values; exp_g_offset <- g_precomp$fitted_values
  reduced_fits <- lapply(seq(1L, n_em_rep), function(i) {
    run_reduced_em_algo(m = m, g = g, exp_m_offset = exp_m_offset, exp_g_offset = exp_g_offset,
                        m_pert_guess = guesses$m_perturbation[i],
                        g_pert_guess = guesses$g_perturbation[i],
                        pi_guess = guesses$pi[i], m_fam = m_fam, g_fam = g_fam)
    })

  # obtain the best run according to log-likelihood
  best_run <- select_best_em_run(reduced_fits)

  # Finally, run full GLM-EIV using pilot estimates
  fit <- run_full_glmeiv_given_pilot_params(m = m, g = g, m_fam = m_fam, g_fam = g_fam,
                                            pi_guess = best_run$pi, m_intercept_guess = m_precomp$fitted_intercept,
                                            m_perturbation_guess = best_run$m_perturbation, m_covariate_coefs_guess = m_precomp$covariate_coefs,
                                            g_intercept_guess = g_precomp$fitted_intercept, g_perturbation_guess = best_run$g_perturbation,
                                            g_covariate_coefs_guess = g_precomp$covariate_coefs, covariate_matrix = covariate_matrix,
                                            m_offset = m_offset, g_offset = g_offset)
  return(fit)
}


#' Run GLM-EIV precomputation
#'
#' Runs a precomputation on a gene or gRNA.
#'
#' @param y a vector of gene or gRNA expressions
#' @param covariate_matrix the matrix of technical factors
#' @param offset a vector of offsets, generally log library size
#' @param fam_str a string specifying the family, currently either "Negative Binomial" or "poisson"
#'
#' @return a list containing (i) estimated theta (in NB case only), (ii) fitted values, (iii) fitted intercept term, (iv) covariate_coef terms
#' @export
run_glmeiv_precomputation <- function(y, covariate_matrix, offset, fam_str) {
  out <- list()
  if (fam_str == "poisson") {
    fit_precomp <- stats::glm(formula = y ~ ., family = poisson(), data = dplyr::mutate(data.frame(y = y), covariate_matrix), offset = offset)
  } else if (fam_str == "Negative Binomial") {
    form <- stats::as.formula(if (is.null(offset)) "y ~ ." else "y ~ . + offset(offset)")
    fit_precomp <- MASS::glm.nb(formula = form, data = dplyr::mutate(data.frame(y = y), covariate_matrix))
    out$theta <- fit_precomp$theta
  } else {
    stop("Family string not recognized.")
  }
  coefs <- as.list(stats::coef(fit_precomp))
  fitted_values <- as.numeric(stats::fitted.values(fit_precomp))
  fitted_intercept <- coefs[["(Intercept)"]]
  covarite_coefs <- unlist(coefs[names(coefs) != "(Intercept)"])
  out$fitted_values <- fitted_values
  out$fitted_intercept <- fitted_intercept
  out$covariate_coefs <- covarite_coefs
  return(out)
}
