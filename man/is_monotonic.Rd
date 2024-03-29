% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{is_monotonic}
\alias{is_monotonic}
\title{check em runs}
\usage{
is_monotonic(em_runs)
}
\arguments{
\item{em_runs}{a list of em runs}
}
\value{
TRUE or FALSE
}
\description{
Verify that the log-likelihood was monotonically increasing across every EM run.
}
