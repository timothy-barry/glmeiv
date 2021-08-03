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


#' Run thresholding method simulatr V2
#'
#' Runs the thresholding method
#'
#' @inheritParams run_full_glmeiv_given_weights
#' @return coefficient table
#' @export
run_thresholding_method_simulatr <- function(dat, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha) {
  # pull g_fam and m_fam from dat, if available
  if ("m_precomp" %in% names(attributes(dat)) && "g_precomp" %in% names(attributes(dat))) {
    m_precomp <- attr(dat, "m_precomp"); m_fam <- m_precomp$fam
    g_precomp <- attr(dat, "g_precomp"); g_fam <- g_precomp$fam
  }
  m <- dat$m
  g <- dat$g
  n <- length(m)
  # run method with timing
  time <- system.time({
    # get the Bayes-optimal threshold
    bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
    # threshold the counts
    phat <- as.integer(g > bdy)
    # check if OK
    s_phat <- sum(phat)
    lower_try_thresh <- (n * pi)/20
    upper_try_thresh <- n - (n * pi)/20
    if (s_phat <= lower_try_thresh || s_phat >= upper_try_thresh) { # too unbalanced; do not attempt fit
      out <- data.frame(parameter = "meta",
                        target = c("fit_attempted", "sum_phat"),
                        value = c(0, s_phat))
    } else {
      # obtain the data matrix
      data_mat <- data.frame(m = m, perturbation = phat) %>% dplyr::mutate(covariate_matrix)
      # fit the model
      fit <- stats::glm(formula = m ~ ., data = data_mat, family = m_fam, offset = m_offset)
      # get the summary
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
        dplyr::add_row(parameter = "meta", target = "sum_phat", value = s_phat)
    }
  })[["elapsed"]]
  # append the run time to out
  out <- out %>% dplyr::add_row(.,parameter = "meta", target = "time", value = time)
  return(out)
}
