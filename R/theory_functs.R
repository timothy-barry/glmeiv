#' Get bias of thresholding estimator
#'
#' Assumes without loss of generality that m_intercept and g_intercept are 0. Also, assumes the Bayes-optimal threshold has been chosen.
#'
#' Assumes errors are N(0,1).
#'
#' @param m_perturbation perturbation coefficient for mRNA model
#' @param g_perturbation perturbation coefficient for gRNA model
#' @param pi probability of perturbation
#'
#' @return the bias of the thresholding estimator
#' @export
get_tresholding_estimator_bias <- function(m_perturbation, g_perturbation, pi) {
  c <- (1/g_perturbation) * (log(1-pi) - log(pi)) + (1/2) * g_perturbation
  z <- 1 - pnorm(c)
  w <- 1 - pnorm(c - g_perturbation)
  mean_phat <- z * (1 - pi) + w * pi
  mult_factor <- exp(log(pi) + log(w - mean_phat) - log(mean_phat) - log(1 - mean_phat))
  bias <- m_perturbation * (1 - mult_factor) * -1
  return(bias)
}
