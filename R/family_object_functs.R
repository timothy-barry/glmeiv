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
  fam_str <- gsub(pattern = '\\([0-9]*\\)', replacement = "", x = f$family)
  new_f <- switch(EXPR = fam_str,
         "poisson" =  augment_poisson_family_object(f),
         "Negative Binomial" = augment_negbinom_family_object(f),
         "gaussian" = augment_gaussian_family_object(f))
  if (is.null(new_f)) stop("Family not recognized. Please use one of the programmed families.")
  new_f$augmented <- TRUE
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
  if (f$link == "log") {
    theta <- get(x = ".Theta", envir = environment(f$variance))
    f$skewness <- function(mu) (2 * mu + theta)/(sqrt(theta * mu) * sqrt(mu + theta))
    f$mu.eta.prime <- function(eta) pmax(exp(eta), .Machine$double.eps)
    f$simulate_from_mus <- function(mus) sapply(X = mus, FUN = function(mu) MASS::rnegbin(n = 1, mu = mu, theta = theta))
    f$log_py_given_mu <- function(y, mu) y * log(mu) - y * log(theta + mu) - theta * log(mu + theta)
  } else stop(get_nonstandard_link_msg())
  f$bayes_classifier <- function(mu0, mu1, pi) (theta * (log(mu0 + theta) - log(mu1 + theta)) + log(pi) - log(1 - pi))/(log(mu0 * (mu1 + theta)) - log(mu1 * (mu0 + theta)))
  f$get_log_lik <- function(object) glm_log_lik(object)
  return(f)
}


augment_poisson_family_object <- function(f) {
  if (f$link == "log") {
    f$skewness <- function(mu) mu^(-1/2)
    f$mu.eta.prime <- function(eta) pmax(exp(eta), .Machine$double.eps)
    f$simulate_from_mus <- function(mus) sapply(X = mus, FUN = function(mu) stats::rpois(1, mu))
    f$log_py_given_mu <- function(y, mu) stats::dpois(x = y, lambda = mu, log = TRUE)
  } else stop(get_nonstandard_link_msg())
  f$bayes_classifier <- function(mu0, mu1, pi) (mu0 - mu1 + log(pi) - log(1 - pi))/(log(mu0) - log(mu1))
  f$get_log_lik <- function(object) glm_log_lik(object)
  return(f)
}


augment_gaussian_family_object <- function(f) {
  if (f$link == "identity") {
    f$skewness <- function(mu) rep(0, length(mu))
    f$mu.eta.prime <- function(eta) rep(0, length(eta))
    f$simulate_from_mus <- function(mus) sapply(X = mus, FUN = function(mu) stats::rnorm(1, mu))
    f$log_py_given_mu <- function(y, mu) -(1/2) * (y - mu)^2
  } else stop(get_nonstandard_link_msg())
  f$bayes_classifier <- function(mu0, mu1, pi) ((1/2) * (mu0^2 - mu1^2) + log(pi) - log(1 - pi))/(mu0 - mu1)
  f$get_log_lik <- function(object) lm_log_lik(object)
  return(f)
}
