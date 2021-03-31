#' Check column names
#'
#' Checks if "perturbation" is a column name of covariate_matrix; if so, produces an error.
#'
#' @param covariate_matrix a matrix of observed covariates
#'
#' @return NULL
check_col_names <- function(covariate_matrix) {
  if ("perturbation" %in% colnames(covariate_matrix)) {
    stop("The column name \"perturbation\" is reserved. Please choose another column name.")
  }
}


#' Diagonal matrix multiply
#'
#' Performs the matrix multiplication ADB, where D is diagonal (and A and B are
#' not necessarily diagonal).
#'
#' @param A a matrix
#' @param D a diagonal matrix, represented as a numeric vector
#' @param B a matrix
#'
#' @return the matrix product ADB.
#' @examples
#' \dontrun{
#' B <- matrix(sample(0:1, size = 12, TRUE), ncol = 3)
#' D <- c(2, 3, 4, 5)
#' A <- t(B)
#' all.equal(A %*% diag(D) %*% B, diag_matrix_multiply(A, D, B))
#' }
diag_matrix_multiply <- function(A, D, B) {
  A %*% (B * D)
}
