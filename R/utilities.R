#' Check column names
#'
#' Checks if "perturbation" is a column name of covariate_matrix; if so, produces an error.
#'
#' @param covariate_matrix
#'
#' @return NULL
#' @noRd
check_col_names <- function(covariate_matrix) {
  if ("perturbation" %in% colnames(covariate_matrix)) {
    stop("The column name \"perturbation\" is reserved. Please choose another column name.")
  }
}
