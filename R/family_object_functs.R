#' Augment family object
#'
#' Augments a family object, adding the fields "skewness," "mu.eta.prime, "simulate_from_mus," and "log_py_given_mu."
#'
#' @param f a family object
#'
#' @return the augmented family object
#' @export
#' @examples
#' pois_augmented <- augment_family_object(poisson())
augment_family_object <- function(f) {
  fam_str <- gsub(pattern = '\\(.*)', replacement = "", x = f$family)
  new_f <- switch(EXPR = fam_str,
         "poisson" =  augment_poisson_family_object(f),
         "Negative Binomial" = augment_negbinom_family_object(f),
         "gaussian" = augment_gaussian_family_object(f))
  if (is.null(new_f)) stop("Family not recognized. Please use one of the programmed families.")
  new_f$augmented <- TRUE
  new_f$fam_str <- fam_str
  return(new_f)
}


#' get_nonstandard_link_msg
#'
#' Produce error message for nonstandard link function.
#'
#' @return NULL
#' @noRd
get_nonstandard_link_msg <- function() {
  "You are using a non-standard link function. Please extend family object manually or use a standard link."
}


#' Augment family object helpers
#'
#' These helper functions augment the family object for Poisson, NB, and Gaussian families.
#'
#' @name augment_family_object_helpers
#' @return an augmented family object
#' @noRd
NULL


augment_negbinom_family_object <- function(f) {
  theta <- get(x = ".Theta", envir = environment(f$variance))
  f$skewness <- function(mu) (2 * mu + theta)/(sqrt(theta * mu) * sqrt(mu + theta))
  f$mu.eta.prime <- function(eta) pmax(exp(eta), .Machine$double.eps)
  f$simulate_from_mus <- function(mus) sapply(X = mus, FUN = function(mu) MASS::rnegbin(n = 1, mu = mu, theta = theta))
  f$simulate_n_times_given_mu <- function(n, mu) MASS::rnegbin(n = n, mu = mu, theta = theta)
  f$log_py_given_mu <- function(y, mu) stats::dnbinom(x = y, size = theta, mu = mu, log = TRUE)
  f$bayes_classifier <- function(mu0, mu1, pi) (theta * (log(mu0 + theta) - log(mu1 + theta)) + log(pi) - log(1 - pi))/(log(mu0 * (mu1 + theta)) - log(mu1 * (mu0 + theta)))
  f$flexmix_fam <- "poisson"
  f$density <- function(mu, xgrid) stats::dnbinom(x = xgrid, size = theta, mu = mu)
  f$d_log_py <- function(y, mu_0, mu_1) y * log(mu_0) - y * log(theta + mu_0) - theta * log(mu_0 + theta) - (y * log(mu_1) - y * log(theta + mu_1) - theta * log(mu_1 + theta))
  f$weighted_log_lik <- function(y, mu_0, mu_1, Ti1s) sum((1 - Ti1s) * stats::dnbinom(x = y, size = theta, mu = mu_0, log = TRUE) + Ti1s * stats::dnbinom(x = y, size = theta, mu = mu_1, log = TRUE))
  return(f)
}


augment_poisson_family_object <- function(f) {
  f$skewness <- function(mu) mu^(-1/2)
  f$mu.eta.prime <- function(eta) pmax(exp(eta), .Machine$double.eps)
  f$simulate_from_mus <- function(mus) sapply(X = mus, FUN = function(mu) stats::rpois(1, mu))
  f$simulate_n_times_given_mu <- function(n, mu) stats::rpois(n, mu)
  f$log_py_given_mu <- function(y, mu) stats::dpois(x = y, lambda = mu, log = TRUE)
  f$bayes_classifier <- function(mu0, mu1, pi) (mu0 - mu1 + log(pi) - log(1 - pi))/(log(mu0) - log(mu1))
  f$flexmix_fam <- "poisson"
  f$density <- function(mu, xgrid) stats::dpois(x = xgrid, lambda = mu)
  f$d_log_py <- function(y, mu_0, mu_1) y * (log(mu_0) - log(mu_1)) + mu_1 - mu_0
  f$weighted_log_lik <- function(y, mu_0, mu_1, Ti1s) sum((1 - Ti1s) * stats::dpois(y, mu_0, TRUE) + Ti1s * stats::dpois(y, mu_1, TRUE))
  return(f)
}


augment_gaussian_family_object <- function(f) {
  f$skewness <- function(mu) rep(0, length(mu))
  f$mu.eta.prime <- function(eta) rep(0, length(eta))
  f$simulate_from_mus <- function(mus) sapply(X = mus, FUN = function(mu) stats::rnorm(1, mu))
  f$simulate_n_times_given_mu <- function(n, mu) stats::rnorm(n, mu)
  f$log_py_given_mu <- function(y, mu) stats::dnorm(x = y, mean = mu, log = TRUE)
  f$bayes_classifier <- function(mu0, mu1, pi) ((1/2) * (mu0^2 - mu1^2) + log(pi) - log(1 - pi))/(mu0 - mu1)
  f$get_log_lik <- function(object) lm_log_lik(object)
  f$flexmix_fam <- "gaussian"
  f$density <- function(mu, xgrid) stats::dnorm(x = xgrid, mean = mu)
  f$d_log_py <- function(y, mu_0, mu_1) -(1/2) * ((y - mu_0)^2 - (y - mu_1)^2)
  f$weighted_log_lik <- function(y, mu_0, mu_1, Ti1s) sum((1 - Ti1s) * stats::dnorm(y, mu_0, log = TRUE) + Ti1s * stats::dpois(y, mu_1, log = TRUE))
  return(f)
}
