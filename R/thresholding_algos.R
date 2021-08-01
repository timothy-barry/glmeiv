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
#' \dontrun{
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 1
#' pi <- 0.2; n <- 1000; B <- 5; alpha <- 0.95
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
#' }
run_thresholding_method_simulatr <- function(dat_list, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha) {
  # first, obtain the optimal boundary
  bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
  n_datasets <- length(dat_list)
  n <- nrow(dat_list[[1]])
  lower_try_thresh <- (pi * n/20); upper_try_thresh <- n - (pi * n/20)
  res_list <- lapply(X = seq(1, n_datasets), FUN = function(i) {
    dat <- dat_list[[i]]
    g <- dat$g
    phat <- as.integer(g > bdy)
    s_phat <- sum(phat)
    if (s_phat <= lower_try_thresh || s_phat >= upper_try_thresh) { # too unbalanced; do not attempt fit
      out <- data.frame(parameter = "meta",
                        target = c("fit_attempted", "sum_phat"),
                        value = c(0, s_phat),
                        run_id = i)
    } else {
      # next, create the data matrix
      data_mat <- data.frame(m = dat$m, perturbation = phat) %>% dplyr::mutate(covariate_matrix)
      # fit model
      fit <- stats::glm(formula = m ~ ., data = data_mat, family = m_fam, offset = m_offset)
      # get the effect size estimates and standard errors
      s <- summary(fit)$coefficients
      row.names(s)[row.names(s) == "(Intercept)"] <- "intercept"
      mult_factor <- stats::qnorm(1 - (1 - alpha)/2)
      confint_lower <- s[,"Estimate"] - mult_factor * s[,"Std. Error"]
      confint_higher <- s[,"Estimate"] + mult_factor * s[,"Std. Error"]
      names(confint_lower) <- names(confint_higher) <- NULL
      out <- data.frame(parameter = paste0("m_", row.names(s)),
                 estimate = s[,"Estimate"],
                 std_error = s[,"Std. Error"],
                 p_value = if (m_fam$family == "poisson") s[,"Pr(>|z|)"] else s[,"Pr(>|t|)"],
                 confint_lower = confint_lower,
                 confint_higher = confint_higher) %>%
        tidyr::pivot_longer(cols = -parameter, names_to = "target") %>%
        dplyr::add_row(parameter = "meta", target = "fit_attempted", value = 1) %>%
        dplyr::add_row(parameter = "meta", target = "sum_phat", value = s_phat) %>%
        dplyr::mutate(run_id = i)
    }
    return(out)
  })
  return(do.call(rbind, res_list))
}


#' Get optimal threshold
#'
#' Thresholds counts according to the Bayes-optimal decision boundary. When a covariate matrix is present, the decision boundary that is used is the mean of the Bayes-optimal decision boundaries across all examples.
#'
#' @param g_intercept intercept for gRNA model
#' @param g_perturbation perturbation coefficient for gRNA model
#' @param g_fam family object describing distribution of g
#' @param pi probability of perturbation
#' @param covariate_matrix (optional) matrix of technical covariates
#' @param g_covariate_coefs (optional) coefficients corresponding to technical factors
#' @param g_offset (optional) the offset vector
#'
#' @return an integer (0/1) vector of thresholded counts.
#' @export
#'
#' @examples
#' \dontrun{
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 1
#' pi <- 0.2; n <- 1000; B <- 500; alpha <- 0.95
#' m_offset <- g_offset <- NULL
#' bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi)
#' }
get_optimal_threshold <- function(g_intercept, g_perturbation, g_fam, pi, covariate_matrix = NULL, g_covariate_coefs = NULL, g_offset = NULL) {
  if (is.null(g_fam$augmented)) g_fam <- augment_family_object(g_fam)
  # get theoretical conditional means
  conditional_means <- compute_theoretical_conditional_means(g_intercept, g_perturbation, g_fam, covariate_matrix, g_covariate_coefs, g_offset)
  # compute the Bayes-optimal boundary
  if (is.null(covariate_matrix)) {
    bdy <- g_fam$bayes_classifier(conditional_means$mu0, conditional_means$mu1, pi)
  } else {
    mu0s <- conditional_means$mu0
    mu1s <- conditional_means$mu1
    bdy <- mean(g_fam$bayes_classifier(mu0s, mu1s, pi))
  }
  return(bdy)
}

