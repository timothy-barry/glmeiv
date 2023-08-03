#' Run two step method
#'
#' Runs the method in which we (1) fit a mixture model to the grna modality, (2) assign perturbation identities based on posterior inclusion probabilities, and then (3) run the GLM
#'
#' @inheritParams run_full_glmeiv_given_weights
#' @return coefficient table
#' @export
run_two_step_method_simulatr <- function(dat, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha) {
  if (!is(m_fam, "family")) m_fam <- m_fam[[1]]
  if (!is(g_fam, "family")) g_fam <- g_fam[[1]]
  # pull g_fam and m_fam from dat, if available
  if ("m_precomp" %in% names(attributes(dat)) && "g_precomp" %in% names(attributes(dat))) {
    m_precomp <- attr(dat, "m_precomp"); m_fam <- m_precomp$fam
    g_precomp <- attr(dat, "g_precomp"); g_fam <- g_precomp$fam
  }
  m <- dat$m
  g <- dat$g
  n <- length(m)

  # fit a mixture model of the grna counts
  g_fit <- flexmix::flexmix(formula = formula(g ~ offset(g_offset)), k = 2, model = flexmix::FLXMRglm(family = "poisson"))
}

