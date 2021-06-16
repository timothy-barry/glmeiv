#' Flip weights
#'
#' @param w a binary (0/1) vector
#' @param p the expected fraction of weights to flip
#'
#' @return a new binary vector with E(p) weights flipped.
#'
#' @examples
#' w <- rbinom(n = 100, size = 1, prob = 0.5)
#' p <- 0.1
#' glmeiv:::flip_weights(w, p)
flip_weights <- function(w, p) {
  out <- (w + stats::rbinom(n = length(w), size = 1, prob = p)) %% 2
  return(out)
}

#' Get optimal threshold
#'
#' Thresholds counts according to the Bayes-optimal decision boundary. When a covariate matrix is present, the decision boundary that is used is the mean of the Bayes-optimal decision boundaries across all examples.
#'
#' @param g_intercept intercept for gRNA model
#' @param g_perturbation perturbation coefficient for gRNA model
#' @param g_fam family object describing distribution of g
#' @param pi probability of perturbation
#' @param covariate_matrix (optional) matrix of technical covariates
#' @param g_covariate_coefs (optional) coefficients corresponding to technical factors
#' @param g_offset (optional) the offset vector
#'
#' @return an integer (0/1) vector of thresholded counts.
#' @export
#'
#' @examples
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 1
#' pi <- 0.2; n <- 1000; B <- 500; alpha <- 0.95
#' m_offset <- g_offset <- NULL
#' bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi)
get_optimal_threshold <- function(g_intercept, g_perturbation, g_fam, pi, covariate_matrix = NULL, g_covariate_coefs = NULL, g_offset = NULL) {
  if (is.null(g_fam$augmented)) g_fam <- augment_family_object(g_fam)
  # get theoretical conditional means
  conditional_means <- compute_theoretical_conditional_means(g_intercept, g_perturbation, g_fam, covariate_matrix, g_covariate_coefs, g_offset)
  # compute the Bayes-optimal boundary
  if (is.null(covariate_matrix)) {
    bdy <- g_fam$bayes_classifier(conditional_means$mu0, conditional_means$mu1, pi)
  } else {
    mu0s <- conditional_means$mu0
    mu1s <- conditional_means$mu1
    bdy <- mean(sapply(seq(1, nrow(covariate_matrix)), function(i) {
      g_fam$bayes_classifier(mu0s[i], mu1s[i], pi)
    }))
  }
  return(bdy)
}


#' Random initialization
#'
#' Output a vector of length n; n * pi of the entries are randomly set to 1, all others to 0.
#'
#' @param n length of output vector
#' @param pi fraction of entries to set to 1
#'
#' @return randomly initialized vector
#' @export
random_initialization <- function(n, pi) {
  out <- integer(n)
  out[sample(x = seq(1, n), size = floor(pi * n), replace = FALSE)] <- 1L
  return(out)
}
