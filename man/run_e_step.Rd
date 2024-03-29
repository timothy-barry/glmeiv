% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/full_em_algo_functs.R
\name{run_e_step}
\alias{run_e_step}
\title{Run E step}
\usage{
run_e_step(
  m_fam,
  g_fam,
  m,
  g,
  m_mus_pert0,
  m_mus_pert1,
  g_mus_pert0,
  g_mus_pert1,
  fit_pi,
  use_mrna_modality = TRUE
)
}
\arguments{
\item{m_fam}{family object describing mRNA counts}

\item{g_fam}{family object describing gRNA counts}

\item{m}{mRNA counts}

\item{g}{gRNA counts}

\item{m_mus_pert0}{fitted (or hypothesized) means mu^m(0)}

\item{m_mus_pert1}{fitted (or hypothesized) means mu^m(1)}

\item{g_mus_pert0}{fitted (or hypothesized) means mu^g(0)}

\item{g_mus_pert1}{fitted (or hypothesized) means mu^g(1)}

\item{fit_pi}{fitted (or hypothesized) value for pi}
}
\value{
a list containing (i) the membership probabilities (Ti1s), (ii) the model log-likelihood, and (iii) the new value of pi (computed ahead of the subsequent M step for convenience).
}
\description{
Runs the E step of the GLM-EIV EM algorithm. Also, returns the updated value of pi (thereby running a portion of the subsequent M step).
}
\details{
Ensures that the value of pi in the next M iteration will be less than 0.5. This core function is called in all variants of GLM-EIV (weights/estimates first, full/reduced, etc.).
}
