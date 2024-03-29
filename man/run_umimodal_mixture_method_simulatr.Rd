% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unimodal_mixture_model_functs.R
\name{run_umimodal_mixture_method_simulatr}
\alias{run_umimodal_mixture_method_simulatr}
\title{Run unimodal mixture method in simulatr}
\usage{
run_umimodal_mixture_method_simulatr(
  dat,
  g_intercept,
  g_perturbation,
  g_fam,
  m_fam,
  pi,
  covariate_matrix,
  g_covariate_coefs,
  m_offset,
  g_offset,
  alpha = 0.95,
  n_em_rep = 15,
  pi_guess_range = c(1e-05, 0.03),
  g_perturbation_guess_range = log(c(0.5, 10)),
  rm_covariate = ""
)
}
\description{
Run unimodal mixture method in simulatr
}
\examples{
m_fam <- g_fam <- augment_family_object(poisson())
n <- 50000
lib_size <- rpois(n = n, lambda = 5000)
m_offset <- g_offset <- log(lib_size)
pi <- 0.01
m_intercept <- log(0.01)
m_perturbation <- log(0.5)
g_intercept <- log(0.005)
g_perturbation <- log(2.5)
covariate_matrix <- data.frame(batch = rbinom(n = n, size = 1, prob = 0.5))
m_covariate_coefs <- log(0.9)
g_covariate_coefs <- log(1.1)
dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept,
m_perturbation = m_perturbation, g_fam = g_fam, g_intercept = g_intercept,
g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
covariate_matrix = covariate_matrix, m_covariate_coefs = m_covariate_coefs,
g_covariate_coefs = g_covariate_coefs, m_offset = m_offset, g_offset = g_offset)[[1]]
m <- dat$m; g <- dat$g; p <- dat$p
fit <- run_umimodal_mixture_method_simulatr(dat, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset)
}
