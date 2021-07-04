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
#' \dontrun{
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 1
#' pi <- 0.2; n <- 1000; B <- 500; alpha <- 0.95
#' m_offset <- g_offset <- NULL
#' bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi)
#' }
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


#' Initialize weights using marginal mixtures
#'
#' Initializes weights for GLM-EIV by fitting marginal mixtures, then taking weighted average.
#'
#' @param m mRNA counts
#' @param g gRNA counts
#' @param m_fam family describing mRNA counts
#' @param g_fam family describing gRNA counts
#' @param m_offset optional offsets for m
#' @param g_offset optional offsets for g
#' @param lambda optional weight in weighted average; if not supplied, chosen adaptively
#'
#' @return initial weights for algorithm
initialize_weights_using_marginal_mixtures <- function(m, g, m_fam, g_fam, m_offset, g_offset, lambda = NULL) {
  m_weights <- get_marginal_mixture_weights(m, m_fam, m_offset)
  g_weights <- get_marginal_mixture_weights(g, g_fam, g_offset)
  if (is.null(lambda)) {
    d_m <- compute_mean_distance_from_half(m_weights)
    d_g <- compute_mean_distance_from_half(g_weights)
    d_sum <- d_m + d_g
    lambda <- if (d_sum < 1e-10) 1 else d_g/(d_sum)
  }
  out <- lambda * g_weights + (1 - lambda) * m_weights
  return(out)
}


#' Get marginal mixture weights
#'
#' Fit a marginal mixture model; get the (correctly oriented) weights
#'
#' @param v a vector (of mRNA or gRNA counts)
#' @param fam family object describing the counts
#' @param offset optional vector of linear offsets
#'
#' @return the EM weights
get_marginal_mixture_weights <- function(v, fam, offset) {
  flex_fit <- flexmix::flexmix(v ~ 1, k = 2,
                               model = flexmix::FLXglm(family = fam$flexmix_fam,
                                                       offset = offset))
  if (flex_fit@k == 1) {
    w <- rep(0.5, length(v))
  } else {
    w_matrix <- flex_fit@posterior$scaled
    w <- if (sum(w_matrix[,1]) <= sum(w_matrix[,2])) w_matrix[,1] else w_matrix[,2]
  }
  return(w)
}


#' Append noise to weights
#'
#' Adds Gaussian noise to an initial weight vector.
#'
#' @param w initial weight vector
#' @param n_rep number of noisy columns to append to w
#' @param sd standard deviation of noise
#'
#' @return initial weight matrix
append_noise_to_weights <- function(w, n_rep, sd) {
  initial_Ti1_matrix <- replicate(n = n_rep, {
    noise <- stats::rnorm(n = length(w), mean = 0, sd = sd)
    out <- w + noise
    out[out > 1] <- 1; out[out < 0] <- 0
    out
  }) %>% cbind(w, .)
  return(initial_Ti1_matrix)
}


compute_mean_distance_from_half <- function(v) {
  n <- length(v)
  sum((v - 0.5)^2)/n
}
