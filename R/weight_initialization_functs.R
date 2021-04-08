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


