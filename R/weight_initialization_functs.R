#' Flip weights
#'
#' @param w a binary (0/1) vector
#' @param p the expected fraction of weights to flip
#'
#' @return a new binary vector with E(p) weights flipped.
#'
#' @examples
#' \dontrun{
#' w <- rbinom(n = 100, size = 1, prob = 0.5)
#' p <- 0.1
#' flip_weights(w, p)
#' }
flip_weights <- function(w, p) {
  out <- (w + stats::rbinom(n = length(w), size = 1, prob = p)) %% 2
  return(out)
}


#' Threshold counts
#'
#' Thresholds counts according to the Bayes-optimal decision boundary. When a covariate matrix is present, the decision boundary that is used is the mean of the Bayes-optimal decision boundaries across all examples.
#'
#' @param g gRNA counts
#' @param g_intercept intercept for gRNA model
#' @param g_pert perturbation coefficient for gRNA model
#' @param g_fam family object describing distribution of g
#' @param pi probability of perturbation
#' @param covariate_matrix (optional) matrix of technical covariates
#' @param g_covariate_coefs (optional) coefficients corresponding to technical factors
#' @param g_offset (optional) the offset vector
#'
#' @return
#' @export
#'
#' @examples
#' g_intercept <- 2
#' g_pert <- 1
#' g_fam <- poisson()
#' pi <- 0.2
#' covariate_matrix <- NULL
#' g_covariate_coefs <- NULL
#' g_offset <- NULL
#' g <- rpois(1000, 5)
#' threshold_counts(g, g_intercept, g_pert, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
#' covariate_matrix <- data.frame(lib_size = runif(1000))
#' g_covariate_coefs <- 2
#' g_offset <- rpois(1000, 4)
#' threshold_counts(g_intercept, g_pert, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
threshold_counts_optimal <- function(g, g_intercept, g_pert, g_fam, pi, covariate_matrix = NULL, g_covariate_coefs = NULL, g_offset = NULL) {
  if (is.null(g_fam$augmented)) g_fam <- augment_family_object(g_fam)
  # get theoretical conditional means
  conditional_means <- compute_theoretical_conditional_means(g_intercept, g_pert, g_fam, covariate_matrix, g_covariate_coefs, g_offset)
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
  return(as.integer(g >= bdy))
}
