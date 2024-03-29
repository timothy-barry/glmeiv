% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulatr_functs.R
\name{generate_full_data}
\alias{generate_full_data}
\title{Generate full data}
\usage{
generate_full_data(
  m_fam,
  m_intercept,
  m_perturbation,
  g_fam,
  g_intercept,
  g_perturbation,
  pi,
  n,
  B,
  covariate_matrix,
  m_covariate_coefs,
  g_covariate_coefs,
  m_offset,
  g_offset,
  run_mrna_unknown_theta_precomputation = FALSE,
  run_grna_unknown_theta_precomputation = FALSE,
  rm_covariate = "",
  mrna_duplet_rate = 0,
  grna_duplet_rate = 0
)
}
\arguments{
\item{m_fam}{family for m}

\item{m_intercept}{m intercept}

\item{m_perturbation}{m perturbation coefficient}

\item{g_fam}{family for g}

\item{g_intercept}{g intercept}

\item{g_perturbation}{g perturbation coefficient}

\item{pi}{probability of perturbation}

\item{n}{number of cells}

\item{B}{number of i.i.d. datasets to generate}

\item{covariate_matrix}{the (fixed) covariate matrix of confounding factors}

\item{m_covariate_coefs}{coefficients for technical factors in model for m}

\item{g_covariate_coefs}{coefficients for technical factors in model for g}

\item{m_offset}{optional (fixed) offset vector for m}

\item{g_offset}{optional (fixed) offset vector for g}

\item{run_unknown_theta_precomputation}{optional; if TRUE, runs the negative binomial (unknown theta) precomputation on all datasets and stores the results ("m_precomp," "g_precomp") as attributes.}
}
\value{
a list of length B of synthetic datasets with columns (p, m, g)
}
\description{
Generates B copies of a full GLM-EIV dataset.
}
\examples{
\dontrun{
library(magrittr)
m_fam <- g_fam <- poisson() \%>\% augment_family_object()
pi <- 0.2; n <- 5000; B <- 5
m_intercept <- log(0.01); m_perturbation <- log(0.25); g_intercept <- log(0.005); g_perturbation <- log(2.5)
m_offset <- log(rpois(n = n, lambda = 10000)); g_offset <- log(rpois(n = n, lambda = 5000))
# no covariates
m_covariate_coefs <- g_covariate_coefs <- covariate_matrix <- NULL
dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
g_intercept, g_perturbation, pi, n, B, covariate_matrix, m_covariate_coefs,
g_covariate_coefs, m_offset, g_offset)
# with covariates
covariate_matrix <- data.frame(p_mito = runif(n, 0, 10))
m_covariate_coefs <- -0.1; g_covariate_coefs <- 0.2
dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
g_intercept, g_perturbation, pi, n, B, covariate_matrix, m_covariate_coefs,
g_covariate_coefs, m_offset, g_offset)
}
}
