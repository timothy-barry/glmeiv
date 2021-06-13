#' Run thresholding method for simulatr
#'
#' Runs the thresholding method; this function is meant to be passed to a simulatr object.
#'
#' @param dat_list a list of synthetic data frames; the columns in each data frame should be "m," "p," and "g."
#' @param g_intercept intercept for gRNA model
#' @param g_perturbation perturbation coefficient for gRNA model
#' @param g_fam family object describing gRNA distribution
#' @param m_fam family object describing mRNA distribution
#' @param pi fraction of cells with perturbation
#' @param covariate_matrix the data frame of cell-specific covariates
#' @param g_covariate_coefs the covariate coefficients for the gRNA model
#' @param m_offset offsets for mRNA model
#' @param g_offset offsets for gRNA model
#' @param alpha confidence intervals are at level (1-alhpa)%
#'
#' @return a data frame with columns parameter, target, value, and run_id
#' @export
#'
#' @examples
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 1
#' pi <- 0.2; n <- 1000; B <- 500; alpha <- 0.95
#' m_offset <- g_offset <- NULL
#' m_covariate_coefs <- g_covariate_coefs <- covariate_matrix <- NULL
#' dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
#' g_intercept, g_perturbation, pi, n, B, covariate_matrix, m_covariate_coefs,
#' g_covariate_coefs, m_offset, g_offset)
#' res <- run_thresholding_method_simulatr(dat_list= dat_list,
#' g_intercept = g_intercept,
#' g_perturbation = g_perturbation,
#' g_fam = g_fam,
#' m_fam = m_fam,
#' pi = pi,
#' covariate_matrix = NULL,
#' g_covariate_coefs = NULL,
#' m_offset = NULL,
#' g_offset = NULL,
#' alpha = alpha)
run_thresholding_method_simulatr <- function(dat_list, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha) {
  # first, obtain the optimal boundary
  bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
  n_datasets <- length(dat_list)
  res_list <- lapply(X = seq(1, n_datasets), FUN = function(i) {
    dat <- dat_list[[i]]
    phat <- as.integer(dat$g > bdy)
    # next, create the data matrix
    data_mat <- data.frame(m = dat$m, perturbation = phat) %>% dplyr::mutate(covariate_matrix)
    # fit model
    fit <- stats::glm(formula = m ~ ., data = data_mat, family = m_fam, offset = m_offset)
    # get the effect size estimates and standard errors
    s <- summary(fit)$coefficients
    row.names(s)[row.names(s) == "(Intercept)"] <- "intercept"
    cis <- suppressMessages(stats::confint(fit, level = alpha))
    data.frame(parameter = paste0("m_", row.names(s)),
               estimate = s[,"Estimate"],
               std_error = s[,"Std. Error"],
               p_value = if (m_fam$family == "gaussian") s[,"Pr(>|t|)"] else s[,"Pr(>|z|)"],
               confint_lower = cis[,1],
               confint_higher = cis[,2]) %>%
      tidyr::pivot_longer(cols = -parameter, names_to = "target") %>%
      dplyr::mutate(run_id = i)
  })
  return(do.call(rbind, res_list))
}
