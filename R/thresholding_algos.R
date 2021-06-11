run_thresholding_method_simulatr <- function(dat, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha) {
  # first, obtain the phats
  phat <- threshold_counts_optimal(dat$g, g_intercept, g_perturbation, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
  # next, create the data matrix
  data_mat <- data.frame(m = dat$m, perturbation = phat) %>% dplyr::mutate(covariate_matrix)
  # fit model
  fit <- stats::glm(formula = m ~ ., data = data_mat, family = m_fam, offset = m_offset)
  # get the effect size estimates and standard errors
  s <- summary(fit)$coefficients
  row.names(s)[row.names(s) == "(Intercept)"] <- "intercept"
  cis <- suppressMessages(stats::confint(fit, level = alpha))
  out <- data.frame(parameter = paste0("m_", row.names(s)),
                    estimate = s[,"Estimate"],
                    std_error = s[,"Std. Error"],
                    p_value = if (m_fam$family == "poisson") s[,"Pr(>|z|)"] else s[,"Pr(>|t|)"],
                    confint_lower = cis[,1],
                    confint_higher = cis[,2]) %>%
  tidyr::pivot_longer(cols = -parameter, names_to = "target")
}
