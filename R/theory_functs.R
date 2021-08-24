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
  z <- 1 - stats::pnorm(c)
  w <- 1 - stats::pnorm(c - g_perturbation)
  mean_phat <- z * (1 - pi) + w * pi
  mult_factor <- exp(log(pi) + log(w - mean_phat) - log(mean_phat) - log(1 - mean_phat))
  bias <- m_perturbation * (1 - mult_factor) * -1
  return(bias)
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
    li0 <- intercept + as.numeric((m %*% c(intercept, covariate_coefs))) + offset
    li1 <- li0 + perturbation_coef
  }
  # compute the (theoretical) conditional means
  mui0 <- fam$linkinv(li0)
  mui1 <- fam$linkinv(li1)
  return(list(mu0 = mui0, mu1 = mui1))
}
