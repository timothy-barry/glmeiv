#' Generate data from model
#'
#' Generates data from the GLM-EIV model. The model treats the (unobserved) binary perturbations as random and all other
#' covariates as fixed.
#'
#' @param m_fam the augmented family object for M.
#' @param g_fam the augmented family object for G.
#' @param m_coef the coefficients of the GLM for M. The vector should be named; two of the names should be "intercept"
#' and "perturbation," and the other names should coincide with the column names of covariate_matrix.
#' @param g_coef the coefficients of the GLM for G. This vector should be named in a similar way to m_coef.
#' @param pi the probability of perturbation.
#' @param covariate_matrix the matrix of observed covariates; assumed to be fixed. Should be of class data.frame.
#' @param m_offset (optional) offsets for M.
#' @param g_offset (optional) offsets for G.
#' @param n (optional, required if covariate_matrix NULL) the number of samples to generate
#'
#' @return a list containing m, g, and p.
#' @export
#'
#' @examples
#' library(magrittr)
#' n <- 10000
#' m_fam <- augment_family_object(poisson())
#' g_fam <- augment_family_object(poisson())
#' m_offset <- 0
#' g_offset <- 0
#' pi <- 0.2
#' covariate_matrix <- data.frame(p_mito = runif(n = n, 0, 10),
#'                                lib_size = rpois(n = n, lambda = 1))
#' m_coef <- c(-1, -2, 1, 0.5)
#' g_coef <- c(-1, 2, 1, 0.5)
#' generated_data <- generate_data_from_model(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix)
#' # verify that we recover ground truth when p is known
#' covariate_matrix_full <- data.frame(perturbation = generated_data$p) %>%
#' dplyr::mutate(covariate_matrix)
#' fit_m <- glm(formula = generated_data$m ~ ., family = m_fam, data = covariate_matrix_full)
#' fit_g <- glm(formula = generated_data$g ~ ., family = g_fam, data = covariate_matrix_full)
generate_data_from_model <- function(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix, m_offset = NULL, g_offset = NULL, n = NULL) {
  # augment family objects, if necessary
  if (is.null(m_fam$augmented)) m_fam <- augment_family_object(m_fam)
  if (is.null(g_fam$augmented)) g_fam <- augment_family_object(g_fam)
  # set offsets to zero, if necessary
  if (is.null(m_offset)) m_offset <- 0; if (is.null(g_offset)) g_offset <- 0
  # verify column names ok
  check_col_names(covariate_matrix)
  # sample unobserved binary covariate p
  if (is.null(covariate_matrix)) {
    if (is.null(n)) stop("covariate_matrix is null. You must supply n.")
  } else {
    n <- nrow(covariate_matrix)
  }
  p <- stats::rbinom(n = n, size = 1, prob = pi)
  # append p to covariate matrix
  if (is.null(covariate_matrix)) {
    covariate_matrix_augmented <- data.frame(perturbation = p)
  } else {
    covariate_matrix_augmented <- dplyr::mutate(covariate_matrix, perturbation = p) %>% dplyr::select(perturbation, everything())
  }

  m <- simulate_glm_data(coefs = m_coef, fam = m_fam, offsets = m_offset, X = covariate_matrix_augmented)
  g <- simulate_glm_data(coefs = g_coef, fam = g_fam, offsets = g_offset, X = covariate_matrix_augmented)
  return(list(m = m, g = g, p = p))
}


#' Simulate GLM data
#'
#' @param coefs the named coefficient vector; one of the names should be "intercept;" the other names should match the column names of X.
#' @param fam the augmented family object
#' @param offsets the offset vector; can be 0.
#' @param X a data frame
#'
#' @return
#' data simulated from the GLM model
#' @noRd
simulate_glm_data <- function(coefs, fam, offsets, X) {
  formula <- paste0("~", paste0(colnames(X), collapse = " + ")) %>% stats::as.formula()
  cov_model <- stats::model.matrix(formula, data = X)

  # compute the linear portion of the model
  l_is <- as.numeric((cov_model %*% coefs)) + offsets
  mu_is <- fam$linkinv(l_is)
  y <- fam$simulate_from_mus(mu_is)
  return(y)
}




#' Get quick simulated data
#'
#' Quickly generates simulated data.
#'
#' @param n number of examples
#'
#' @return a list containing the entries m, g, p, pi, covariate_matrix, m_coef, and g_coef.
#' @export
get_quick_simulated_data <- function(n = 1000) {
  m_fam <- augment_family_object(stats::poisson())
  g_fam <- augment_family_object(stats::poisson())
  m_offset <- 0
  g_offset <- 0
  pi <- 0.2
  covariate_matrix <- data.frame(p_mito = stats::runif(n = n, 0, 10),
                                 lib_size = stats::rpois(n = n, lambda = 1))
  m_coef <- c(-2, 3, 0.75, -0.5)
  g_coef <- c(-2, 3, 1, 0.5)
  out <- generate_data_from_model(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix)
  out$pi <- pi
  out$covariate_matrix <- covariate_matrix
  out$m_coef <- m_coef
  out$g_coef <- g_coef
  out$m_fam <- m_fam
  out$g_fam <- g_fam
  return(out)
}
