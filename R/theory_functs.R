#' Get bias of thresholding estimator
#'
#' Assumes without loss of generality that m_intercept and g_intercept are 0. Also, assumes the Bayes-optimal threshold has been chosen.
#'
#' Assumes errors are N(0,1).
#'
#' @param m_perturbation perturbation coefficient for mRNA model
#' @param g_perturbation perturbation coefficient for gRNA model
#' @param g_intercept intercept term in gRNA model
#' @param pi probability of perturbation
#' @param c the threshold
#' @param return_bias if TRUE, returns (absolute) bias, defined as ; if FALSE, returns attenuation fraction
#'
#' @return the bias of the thresholding estimator
#' @export
#' @examples
#' m_perturbation <- 1
#' g_perturbation <- 1
#' g_intercept <- 0.5
#' pi <- 0.5
#' c <- 1
#' get_tresholding_estimator_bias(m_perturbation, g_perturbation, g_intercept, pi, c)
get_tresholding_estimator_bias <- function(m_perturbation, g_perturbation, g_intercept, pi, c, return_bias = TRUE) {
  w <- pnorm(g_intercept + g_perturbation - c)
  z <- pnorm(g_intercept - c)
  exp_phat <- z * (1 - pi) + w * pi
  gamma <- exp(log(pi) + log(w - exp_phat) - log(exp_phat) - log(1 - exp_phat))
  if (return_bias) {
    out <- m_perturbation * (1 - gamma)
  } else {
    out <- gamma
  }
  return(out)
}


#' Get variance of thresholding estimator (no intercept)
#'
#' Returns variance of thresholding estimator in no-intercept model.
#'
#' @param c selected threshold
#' @param g_beta perturbation coefficient for gRNA model
#' @param pi probability of perturbation
#' @param m_beta perturbation coefficient for mRNA model
#' @param n number of cells
#'
#' @return theoretical variance of estimator
#' @export
get_thresholding_estimator_var_no_int <- function(c, g_beta, pi, m_beta, n) {
  z <- 1 - stats::pnorm(c)
  w <- 1 - stats::pnorm(c - g_beta)
  mean_phat <- z * (1 - pi) + w * pi
  l <- get_thresholding_estimator_est_no_int(c, g_beta, pi, m_beta)
  var <- (m_beta^2 * w * pi + mean_phat - 2 * l * m_beta * w * pi + l^2 * mean_phat)/(mean_phat^2)
  return(var/n)
}


#' Get thresholding estimator estimate (no intercept model)
#'
#' @param c selected threshold
#' @param g_beta perturbation coefficient for gRNA model
#' @param pi probability of perturbation
#' @param m_beta perturbation coefficient for mRNA model
#'
#' @return limit (in probability) of estimator
#' @export
get_thresholding_estimator_est_no_int <- function(c, g_beta, pi, m_beta) {
  z <- 1 - stats::pnorm(c)
  w <- 1 - stats::pnorm(c - g_beta)
  mean_phat <- z * (1 - pi) + w * pi
  est <- (m_beta * w * pi)/(z * (1 - pi) + w * pi)
  return(est)
}


#' Compute theoretical conditional means
#'
#' Computes the conditional means of the GLM-EIV model given the ground truth.
#'
#' @param intercept intercept term
#' @param perturbation_coef coefficient corresponding to perturbation
#' @param fam family object describing response distribution
#' @param covariate_matrix (optional) matrix of technical covariates
#' @param covariate_coefs (optional) coefficients corresponding to technical covariates
#' @param offset (optional) vector of offsets
#'
#' @return the scalar (or vector, if covariate_matrix is supplied) of conditional means
compute_theoretical_conditional_means <- function(intercept, perturbation_coef, fam, covariate_matrix = NULL, covariate_coefs = NULL, offset = NULL) {
  # augment family object and set offset to 0, if necessary
  if (is.null(offset)) offset <- 0
  # compute the (theoretical) conditional linear components
  if (is.null(covariate_matrix)) {
    li0 <- intercept + offset
    li1 <- li0 + perturbation_coef
  } else {
    # col_class <- sapply(colnames(covariate_matrix), function(colname) class(covariate_matrix[[colname]]))
    form_str <- paste0("~", paste0(colnames(covariate_matrix), collapse = " + "))
    m <- stats::model.matrix(stats::as.formula(form_str), covariate_matrix)
    li0 <- as.numeric((m %*% c(intercept, covariate_coefs))) + offset
    li1 <- li0 + perturbation_coef
  }
  # compute the (theoretical) conditional means
  mui0 <- fam$linkinv(li0)
  mui1 <- fam$linkinv(li1)
  return(list(mu0 = mui0, mu1 = mui1))
}
