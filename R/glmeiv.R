utils::globalVariables(c("n", "perturbation", "everything", "variable", "counts", "..count..", ".", "run", "iteration", "log_lik", "parameter"))

#' glmeiv: a glm-oriented variant of the errors-in-variables model.
#'
#' GLM-EIV is GLM-oriented variant of the errors-in-variables model. GLM-EIV models the response as
#' a GLM of an unobserved predictor. A noisy version of the predictor is obeserved instead, which
#' itself is modeled as a GLM of the true predictor. glmeiv implements estimation and inference in the GLM-EIV model.
#'
#' @importFrom magrittr %>%
#' @docType package
#' @name glmeiv
NULL
