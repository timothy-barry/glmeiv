#' Run inference on EM fit
#'
#' @param em_fit a fitted GLM-EIV model, as outputted by run_em_algo_given_init.
#' @param alpha (1-alpha) confidence interval produced
#'
#' @return a data frame with columns variable, estimate, std_error, p_value, confint_lower, confint_higher
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- get_quick_simulated_data(1000)
#' em_fit <- run_em_algo_given_init(m = dat$m, g = dat$g, m_fam = poisson(), g_fam = poisson(),
#' covariate_matrix = dat$covariate_matrix, initial_Ti1s = dat$p, m_offset = NULL, g_offset = NULL)
#' se_table <- run_inference_on_em_fit(em_fit)
#' }
run_inference_on_em_fit <- function(em_fit, alpha = 0.95) {
  # compute the asymptotic variance-covariance matrix
  fit_m <- em_fit$fit_m
  fit_g <- em_fit$fit_g
  pieces <- get_info_mat_pieces(em_fit)
  info_mat <- compute_info_matrix(pieces)
  info_mat_names <- c("pi", paste0("m_", names(stats::coef(fit_m))), paste0("g_", names(stats::coef(fit_g))))
  colnames(info_mat) <- row.names(info_mat) <- info_mat_names
  vcov_mat <- solve(info_mat)
  # construct standard errors and confidence intervals
  estimate <- c(em_fit$fit_pi, stats::coef(fit_m), stats::coef(fit_g))
  names(estimate) <- colnames(vcov_mat)
  diag_entries <- diag(vcov_mat)
  diag_entries[diag_entries < 0] <- NA
  std_error <- sqrt(diag_entries)
  z_value <- estimate/std_error
  p_value <- 2 * stats::pnorm(-abs(z_value))
  mult_factor <- stats::qnorm(1 - (1 - alpha)/2)
  confint_lower <- estimate - mult_factor * std_error
  confint_higher <- estimate + mult_factor * std_error
  coef_table <- cbind(estimate, std_error, p_value, confint_lower, confint_higher)
  vars <- row.names(coef_table)
  coef_table <- as.data.frame(coef_table) %>% dplyr::mutate(variable = vars) %>%
    dplyr::select(variable, everything())
  row.names(coef_table) <- NULL
  coef_table$variable <- gsub(pattern = "\\(Intercept\\)", replacement = "intercept", x = coef_table$variable)
  return(coef_table)
}


#' Get SE pieces
#'
#' Obtain the pieces needed to compute standard errors.
#'
#' @param em_fit a fitted GLM-EIV model, as outputted by run_em_algo_given_init.
#'
#' @return A list of pieces required to calculate the observed information matrix.
#' @examples
#' \dontrun{
#' dat <- get_quick_simulated_data(10000)
#' em_fit <- run_em_algo_given_init(m = dat$m, g = dat$g, m_fam = poisson(), g_fam = poisson(),
#' covariate_matrix = dat$covariate_matrix, initial_Ti1s = dat$p, m_offset = NULL, g_offset = NULL)
#' info_mat_pieces <- get_info_mat_pieces(em_fit)
#' }
get_info_mat_pieces <- function(em_fit) {
  n <- em_fit$n
  idx_0s <- seq(1, n)
  idx_1s <- seq(n + 1, 2 * n)
  fit_m <- em_fit$fit_m
  fit_g <- em_fit$fit_g

  # First, the data matrix
  dat_full <- stats::model.matrix(fit_m)
  X0 <- dat_full[idx_0s,]
  X1 <- dat_full[idx_1s,]

  # Second, the weights
  weights_full <- fit_m$prior.weights
  Ti0s <- as.numeric(weights_full[idx_0s])
  Ti1s <- as.numeric(weights_full[idx_1s])

  # Third, the m's
  m0 <- get_info_mat_pieces_helper(fit_m, idx_0s)
  m1 <- get_info_mat_pieces_helper(fit_m, idx_1s)

  # Finally, the g's
  g0 <- get_info_mat_pieces_helper(fit_g, idx_0s)
  g1 <- get_info_mat_pieces_helper(fit_g, idx_1s)

  out <- list(zero = list(X = X0, Ti = Ti0s, m = m0, g = g0),
              one = list(X = X1, Ti = Ti1s, m = m1, g = g1),
              n = n, pi = em_fit$fit_pi)
  return(out)
}


#' Get matrices from fit
#'
#' Obtains the matrices Delta, V, Delta', and H, and Delta * H from a fitted GLM object
#' given a set of indexes.
#'
#' @param fit a fitted GLM object
#' @param idxs the sample indexes on which to extract the matrices
#'
#' @return a list of diagonal matrices, represented as list of numeric vectors. The diagonal matrices are
#' Delta, V, Delta', H, and Delta * H.
get_info_mat_pieces_helper <- function(fit, idxs) {
  fam <- fit$family
  l_is <- as.numeric(fit$linear.predictors[idxs])
  y_is <- as.numeric(fit$y[idxs])
  mu_is <- fam$linkinv(l_is)
  var_is <- fam$variance(mu_is)
  h_prime_l_is <- fam$mu.eta(l_is)/var_is
  skewness_is <- fam$skewness(mu_is)
  mu_eta_prime_is <- fam$mu.eta.prime(l_is)
  h_dprime_l_is <- (mu_eta_prime_is - var_is^(3/2) * skewness_is * (h_prime_l_is^2))/var_is
  h <- y_is - mu_is
  return(list(Delta = h_prime_l_is, Delta_prime = h_dprime_l_is, V = var_is, H = h, DeltaH = h_prime_l_is * h))
}


#' Compute information submatrix
#'
#' Compute information submatrices 1-6. For matrices 2-3 and 5-6, specify either "m" or "g."
#'
#' @param pieces "pieces list," as outputted by get_info_mat_pieces.
#' @param k either "m" or "g;" required for compute_info_submat_23 and compute_info_submat_56.
#' @name compute_info_submat
#'
#' @return the requested information submatrix
NULL


#' @rdname compute_info_submat
compute_info_submat_1 <- function(pieces) {
  Ti1 <- pieces[["one"]][["Ti"]]
  n <- pieces$n
  pi_hat <- pieces$pi
  out <- (1/pi_hat^2 - 1/(1 - pi_hat)^2) * sum(Ti1) + n/(1 - pi_hat)^2 + (1/(1 - pi_hat) + 1/pi_hat)^2 * (sum(Ti1^2 -Ti1))
  return(out)
}


#' @rdname compute_info_submat
compute_info_submat_23 <- function(pieces, k) {
  ss <- c("zero", "one")
  g <- dplyr::mutate_all(.tbl = expand.grid(ss = ss, ts = ss), as.character)

  m.1 <- lapply(ss, function(s) {
    D <- pieces[[s]][["Ti"]] * (pieces[[s]][[k]][["Delta"]]^2 * pieces[[s]][[k]][["V"]] - pieces[[s]][[k]][["Delta_prime"]] * pieces[[s]][[k]][["H"]])
    diag_matrix_multiply(t(pieces[[s]][["X"]]), D, pieces[[s]][["X"]])
  })

  m.2 <- mapply(function(s, tau) {
    D <- pieces[[s]][["Ti"]] * pieces[[s]][[k]][["DeltaH"]] * pieces[[tau]][["Ti"]] * pieces[[tau]][[k]][["DeltaH"]]
    diag_matrix_multiply(t(pieces[[s]][["X"]]), D, pieces[[tau]][["X"]])
  }, g$ss, g$ts, SIMPLIFY = FALSE)

  m.3 <- lapply(ss, function(s) {
    D <- pieces[[s]][["Ti"]] * pieces[[s]][[k]][["DeltaH"]]^2
    diag_matrix_multiply(t(pieces[[s]][["X"]]), D, pieces[[s]][["X"]]) * -1
  })

  names(m.1) <- names(m.2) <- names(m.3) <- NULL
  to_sum <- c(m.1, m.2, m.3)
  out <- Reduce('+', to_sum)
  return(out)
}


#' @rdname compute_info_submat
compute_info_submat_4 <- function(pieces) {
  ss <- c("zero", "one")
  g <- dplyr::mutate_all(.tbl = expand.grid(ss = ss, ts = ss), as.character)

  m.1 <- mapply(function(s, tau) {
    D <- pieces[[s]][["Ti"]] * pieces[[s]][["g"]][["DeltaH"]] * pieces[[tau]][["Ti"]] * pieces[[tau]][["m"]][["DeltaH"]]
    diag_matrix_multiply(t(pieces[[s]][["X"]]), D, pieces[[tau]][["X"]])
  }, g$ss, g$ts, SIMPLIFY = FALSE)

  m.2 <- lapply(ss, function(s) {
    D <- pieces[[s]][["Ti"]] * pieces[[s]][["g"]][["DeltaH"]] * pieces[[s]][["m"]][["DeltaH"]]
    diag_matrix_multiply(t(pieces[[s]][["X"]]), D, pieces[[s]][["X"]]) * -1
  })

  names(m.1) <- names(m.2) <- NULL
  to_sum <- c(m.1, m.2)
  out <- Reduce('+', to_sum)
  return(out)
}


#' @rdname compute_info_submat
compute_info_submat_56 <- function(pieces, k) {
  pi_hat <- pieces$pi
  u_k <- pieces[["zero"]][["Ti"]] * pieces[["one"]][["Ti"]] * pieces[["zero"]][[k]][["DeltaH"]]
  v_k <- pieces[["zero"]][["Ti"]] * pieces[["one"]][["Ti"]] * pieces[["one"]][[k]][["DeltaH"]]
  output <- (1/pi_hat + 1/(1 - pi_hat)) * (t(pieces[["zero"]][["X"]]) %*% u_k - t(pieces[["one"]][["X"]]) %*% v_k)
  return(output)
}


compute_info_matrix <- function(pieces) {
  submat_1 <- compute_info_submat_1(pieces)
  submat_2 <- compute_info_submat_23(pieces, "m")
  submat_3 <- compute_info_submat_23(pieces, "g")
  submat_4 <- compute_info_submat_4(pieces)
  submat_5 <- compute_info_submat_56(pieces, "g")
  submat_6 <- compute_info_submat_56(pieces, "m")

  info_mat <- rbind(cbind(submat_1, t(submat_6), t(submat_5)),
                    cbind(submat_6, submat_2, t(submat_4)),
                    cbind(submat_5, submat_4, submat_3))
  return(info_mat)
}
