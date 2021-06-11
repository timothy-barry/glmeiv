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


#' check em runs
#'
#' Verify that the log-likelihood was monotonically increasing across every EM run.
#'
#' @param em_runs a list of em runs
#'
#' @return TRUE or FALSE
is_monotonic <- function(em_runs) {
  monotonically_increasing <- sapply(X = em_runs, function(run) {
    all(diff(run$log_liks) >= -0.1)
  })
  return(monotonically_increasing)
}


#' Select best EM run
#'
#' Returns the best run, defined as the run with positive perturbation coefficient for G and
#' greatest log-likelihood.
#'
#' @param em_runs a list of outputs of run_em_algo_given_init
#'
#' @return the best run of the list
#' @export
select_best_em_run <- function(em_runs) {
  # select the run with greatest likelihood
  log_liks <- sapply(em_runs, function(run) run$log_lik)
  return(em_runs[[which.max(log_liks)]])
}
