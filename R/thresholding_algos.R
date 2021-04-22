#' Run thresholding method
#'
#' Runs the thresholding method.
#'
#' @param m a vector of m counts
#' @param m_fam a family object describing m
#' @param covariate_matrix the covariate matrix
#' @param p_hat the "estimated" vector of perturbations obtained through thresholding
#' @param m_offset (optional) offset term for m
#' @param alpha (optional; default .95) returns an alpha-CI for each estimated parameter
#'
#' @return a data frame of estimated coefficients, p-values, and standard errors
#' @export
#'
#' @examples
#' n <- 1000
#' m_fam <- poisson()
#' g_fam <- poisson()
#' m_offset <- g_offset <- NULL
#' pi <- 0.2
#' covariate_matrix <- NULL
#' m_coef <- c(2, -1)
#' g_coef <- c(-1, 2)
#' generated_data <- generate_data_from_model(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix, n = n)
#' m <- generated_data$m
#' g <- generated_data$g
#' p_hat <- threshold_counts_no_covariates(g_coef[1], g_coef[2], g, g_fam, pi)
#' run_thresholding_method(m, m_fam, covariate_matrix, p_hat)
run_thresholding_method <- function(m, m_fam, covariate_matrix, p_hat, m_offset = NULL, alpha = 0.95) {
  covariate_matrix_full <- data.frame(perturbation = as.integer(p_hat)) %>%
    dplyr::mutate(covariate_matrix)
  fit <- stats::glm(formula = m ~ ., family = m_fam,
                    data = covariate_matrix_full, offset = m_offset)
  s <- summary(fit)$coefficients
  cis <- suppressMessages(stats::confint(fit, level = alpha))
  out <- data.frame(variable = paste0("m_", row.names(s)),
             estimate = s[,"Estimate"],
             std_error = s[,"Std. Error"],
             p_value = if (m_fam$family == "poisson") s[,"Pr(>|z|)"] else s[,"Pr(>|t|)"],
             confint_lower = cis[,1],
             confint_higher = cis[,2])
  row.names(out) <- NULL
  return(out)
}
