run_umimodal_mixture_method_simulatr <- function(dat, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha) {
  if (!is(m_fam, "family")) m_fam <- m_fam[[1]]
  if (!is(g_fam, "family")) g_fam <- g_fam[[1]]
  # pull g_fam and m_fam from dat, if available
  if ("m_precomp" %in% names(attributes(dat))) {
    m_precomp <- attr(dat, "m_precomp"); m_fam <- m_precomp$fam
  }
  if ("g_precomp" %in% names(attributes(dat))) {
    g_precomp <- attr(dat, "g_precomp"); g_fam <- g_precomp$fam
  }

  # step 0: update the data frame to remove m and p; add the covariates
  dat <- dat |> dplyr::mutate(covariate_matrix)

  # step 1: fit the unimodal mixture model to the gRNA modality
  if (is.null(covariate_matrix)) { # no covariates: varying intercept model with offsets
    fit <- flexmix::flexmix(g ~ 1, k = 2, data = dat,
                   model = flexmix::FLXMRglmfix(family = "poisson", offset = g_offset))
  } else { # covariates, varying intercept, offsets
    covariates <- colnames(covariate_matrix)
    my_form <- paste0("~ ", paste0(covariates, collapse = "+")) |> as.formula()
    fit <- flexmix::flexmix(g ~ 1, k = 2, data = dat,
                   model = flexmix::FLXMRglmfix(family = "poisson", offset = g_offset, fixed = my_form))
  }

  # step 2: obtain the cluster assignments; ensure 1 (perturbed) is the smaller cluster and 0 (unperturbed) is the larger one
  phat <- clusters(fit) - 1L
  if (mean(phat) >= 0.5) {
    phat <- 1L - phat
  }

  # step 3: call the thresholding method simulatr helper
  m <- dat$m
  out <- thresholding_method_simulatr_helper(phat, m, pi, covariate_matrix, m_fam, m_offset, alpha)
  return(out)
}
