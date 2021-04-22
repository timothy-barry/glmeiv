#' GLM log lik
#'
#' Internal function; return the log-likelihood of a GLM object.
#'
#' @param object a fitted GLM
#'
#' @return its log-likelihood
glm_log_lik <- function(object) {
p <- object$rank
val <- p - object$aic/2
}


#' LM log lik
#'
#' Internal function; return the log-likelihood of an lm.
#'
#' @param object a linear model (or GLM with Gaussian family, identity link)
#'
#' @return its log-likelihood
lm_log_lik <- function(object) {
  res <- object$residuals
  p <- object$rank
  N <- length(res)
  if (is.null(w <- object$weights)) {
    w <- rep.int(1, N)
  }
  else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }
  N0 <- N
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
return(val)
}
