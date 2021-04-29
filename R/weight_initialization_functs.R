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


#' Get Bayes-optimal boundary
#'
#' Given an augmented family object, linear component corresponding to P = 0,
#' linear component corresponding to P = 1, and pi, get the Bayes-optimal
#' decision boundary.
#'
#' @param l_0 l_i(0)
#' @param l_1 l_i(1)
#' @param pi probability of perturbation
#' @param fam augmented family object
#'
#' @return the Bayes-optimal decision boundary
#' @examples
#' fam <- poisson() %>% augment_family_object()
#' l_0 <- -1
#' l_1 <- 2
#' pi <- 0.2
#' get_bayes_optimal_bdy(l_0, l_1, pi, fam)
get_bayes_optimal_bdy <- function(l_0, l_1, pi, fam) {
  mu_0 <- fam$linkinv(l_0)
  mu_1 <- fam$linkinv(l_1)
  fam$bayes_classifier(mu_0, mu_1, pi)
}


#' Threshold counts (no covariates)
#'
#' Threshold the g counts in the absence of covariates.
#'
#' @param g_intercept (known) intercept of g model
#' @param g_pert (known) perturbation coefficient of g model
#' @param g g counts
#' @param g_fam family object describing g
#' @param pi (known) probability of perturbation
#'
#' @return The optimally thresholded counts p_hat
#' @export
threshold_counts_no_covariates <- function(g_intercept, g_pert, g, g_fam, pi) {
  fam <- g_fam %>% augment_family_object()
  l_0 <- g_intercept
  l_1 <- g_intercept + g_pert
  if (l_0 == l_1) { # no separation; random guess
    p_hat <- sample(x = c(0,1), size = length(g), replace = TRUE)
  } else {
    bdy <- get_bayes_optimal_bdy(l_0, l_1, pi, fam)
    p_hat <- as.integer(g >= bdy)
  }
  return(p_hat)
}
