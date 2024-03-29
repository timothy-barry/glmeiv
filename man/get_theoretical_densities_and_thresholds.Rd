% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functs.R
\name{get_theoretical_densities_and_thresholds}
\alias{get_theoretical_densities_and_thresholds}
\title{Get DFs for mixture plots}
\usage{
get_theoretical_densities_and_thresholds(sim_spec, xgrid)
}
\arguments{
\item{sim_spec}{a simulatr specifier object}

\item{xgrid}{grid over which to compute the density}
}
\value{
a list containing DFs to plot and optimal thresholds for both mRNA and gRNA
}
\description{
Returns DFs to plot the mixture distributions.
}
