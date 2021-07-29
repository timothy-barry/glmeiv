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
#'
#' @examples
#' \dontrun{
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
#' }
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


#' Run GLM-EIV with known p
#'
#' @param m m data vector
#' @param g g data vector
#' @param m_fam family object for m
#' @param g_fam family object for g
#' @param covariate_matrix the matrix of observed covariates
#' @param p the vector of p's; these typically are unobserved but will be used as the initial weights here
#' @param m_offset (default NULL) the vector of offsets for m
#' @param g_offset (default NULL) the vector of offsets for g
#' @param n_runs (default 10) number of EM algo runs
#' @param p_flip (default 0.15) expected fraction of initial weights to flip in each run
#' @param ep_tol (detault 0.1) tolerance threshold for EM convergence
#' @param max_it (default 50) maximum number of EM iterations (per run)
#' @param alpha (default 0.95) confidence interval level
#' @param reduced_output (default TRUE) return only the best EM run (as determined by log-likelihood)?
#'
#' @return the best EM run
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- get_quick_simulated_data(5000)
#' em_coefs <- run_glmeiv_known_p(m = dat$m, g = dat$g, m_fam = dat$m_fam, g_fam = dat$m_fam,
#' covariate_matrix = dat$covariate_matrix, p = dat$p, m_offset = NULL, g_offset = NULL)
#' # dat$m_coef; dat$g_coef; dat$pi
#' }
run_glmeiv_known_p <- function(m, g, m_fam, g_fam, covariate_matrix, p, m_offset = NULL, g_offset = NULL, n_runs = 5, p_flip = 0.15, ep_tol = 0.5 * 1e-4, max_it = 50, alpha = 0.95, reduced_output = TRUE) {
  initial_Ti1_matrix <- replicate(n_runs, expr = flip_weights(p, p_flip))
  em_runs <- run_em_algo_multiple_inits(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1_matrix, m_offset, g_offset, ep_tol = ep_tol, max_it = max_it)
  if (reduced_output) {
    best_run <- select_best_em_run(em_runs)
    out <- run_inference_on_em_fit(best_run, alpha)
  } else {
    out <- em_runs
  }
  return(out)
}


#' Random initialization
#'
#' Output a vector of length n; n * pi of the entries are randomly set to 1, all others to 0.
#'
#' @param n length of output vector
#' @param pi fraction of entries to set to 1
#'
#' @return randomly initialized vector
#' @export
random_initialization <- function(n, pi) {
  out <- integer(n)
  out[sample(x = seq(1, n), size = floor(pi * n), replace = FALSE)] <- 1L
  return(out)
}


#' Run em algorithm for simulatr using optimal threshold initialization
#'
#' @param dat_list a list of data frames, each of which has columns m, g, p.
#' @param g_intercept the intercept for the gRNA model.
#' @param g_perturbation the perturbation coefficient for the gRNA model.
#' @param g_fam family object describing gRNA model.
#' @param m_fam family object describing mRNA model.
#' @param pi probability of perturbation
#' @param covariate_matrix data frame of technical factors; can be null
#' @param g_covariate_coefs technical factor coefficients for gRNA model
#' @param m_offset optional offset vector for mRNA model
#' @param g_offset optional offset vector for gRNA model
#' @param alpha (1-alpha)% CIs returned
#' @param n_em_rep number of EM algorithm runs to conduct
#' @param p_flip probability of flipping a given initialization weight.
#'
#' @return a data frame with columns parameter, target, value, and run_id
#' @export
#'
#' @examples
#' \dontrun{
#' library(magrittr)
#' m_fam <- g_fam <- poisson() %>% augment_family_object()
#' m_intercept <- 2; m_perturbation <- -1; g_intercept <- -2; g_perturbation <- 3
#' pi <- 0.2; n <- 1000; B <- 5; alpha <- 0.95; n_em_rep <- 5; p_flip <- 0.01
#' m_offset <- g_offset <- NULL
#' m_covariate_coefs <- g_covariate_coefs <- covariate_matrix <- NULL
#' dat_list <- generate_full_data(m_fam, m_intercept, m_perturbation, g_fam,
#' g_intercept, g_perturbation, pi, n, B, covariate_matrix,
#' m_covariate_coefs, g_covariate_coefs, m_offset, g_offset)
#' run_em_algo_simulatr_optimal_thresh(dat_list, g_intercept, g_perturbation,
#' g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset,
#' alpha, n_em_rep, p_flip)
#' }
run_em_algo_simulatr_optimal_thresh <- function(dat_list, g_intercept, g_perturbation, g_fam, m_fam, pi, covariate_matrix, g_covariate_coefs, m_offset, g_offset, alpha, n_em_rep, p_flip) {
  # first, obtain the optimal boundary
  bdy <- get_optimal_threshold(g_intercept, g_perturbation, g_fam, pi, covariate_matrix, g_covariate_coefs, g_offset)
  n_datasets <- length(dat_list)
  n <- nrow(dat_list[[1]])
  res_list <- lapply(X = seq(1, n_datasets), FUN = function(i) {
    dat <- dat_list[[i]]
    g <- dat$g
    phat <- as.integer(g > bdy)
    if (all(phat == 1) || all(phat == 0)) { # just randomly initialize instead
      initial_Ti1_matrix <- replicate(n = n_em_rep, random_initialization(n, pi))
    } else {
      initial_Ti1_matrix <- cbind(phat, replicate(n_em_rep - 1, flip_weights(phat, p_flip)))
    }
    em_fit <- run_em_algo_multiple_inits(m = dat$m, g = dat$g, m_fam = m_fam, g_fam = g_fam,
                                         covariate_matrix = covariate_matrix, initial_Ti1_matrix = initial_Ti1_matrix,
                                         m_offset = m_offset, g_offset = g_offset, return_best = TRUE)
    s <- run_inference_on_em_fit(em_fit, alpha) %>% dplyr::rename("parameter" = "variable")
    tidyr::pivot_longer(s, cols = -parameter, names_to = "target") %>%
      dplyr::add_row(parameter = "information", target = "converged", value = em_fit$converged) %>%
      dplyr::mutate(run_id = i)
  })
  return(do.call(rbind, res_list))
}


#' Flip weights
#'
#' @param w a binary (0/1) vector
#' @param p the expected fraction of weights to flip
#'
#' @return a new binary vector with E(p) weights flipped.
#'
#' @examples
#' w <- rbinom(n = 100, size = 1, prob = 0.5)
#' p <- 0.1
#' glmeiv:::flip_weights(w, p)
flip_weights <- function(w, p) {
  out <- (w + stats::rbinom(n = length(w), size = 1, prob = p)) %% 2
  return(out)
}


plot_em_runs <- function(em_runs) {
  to_plot <- lapply(seq(1, length(em_runs)), function(i) {
    curr_log_liks <- em_runs[[i]]$log_liks
    data.frame(log_lik = curr_log_liks, iteration = seq(1, length(curr_log_liks))) %>%
      dplyr::mutate(run = i)
  }) %>% do.call(what = rbind, args = .) %>% dplyr::mutate(run = factor(run))
  ggplot2::ggplot(to_plot %>% dplyr::filter(iteration >= 3), ggplot2::aes(x = iteration, y = log_lik, col = run)) +
    ggplot2::geom_line() + ggplot2::theme_bw()
}


#' Plot count distribution
#'
#' Plots the distribution of m or g counts, colored by perturbation status.
#'
#' @param generated_data a list containing p, m, and/r g.
#' @param modality either "mRNA" or "gRNA"
#'
#' @return a ggplot of the histogram
#' @export
#'
#' @examples
#' \dontrun{
#' n <- 10000
#' m_fam <- g_fam <- poisson()
#' pi <- 0.4
#' covariate_matrix <- NULL
#' m_coef <- c(1, -2)
#' g_coef <- c(-2, 3)
#' generated_data <- generate_data_from_model(m_fam, g_fam,
#' m_coef, g_coef, pi, covariate_matrix, n = n)
#' p <- plot_count_distribution(generated_data, "gRNA")
#' plot(p)
#' }
plot_count_distribution <- function(generated_data, modality) {
  df <- data.frame(p = factor(x = generated_data$p,  levels = c(1, 0), labels = c("Perturbation", "No perturbation")),
                   counts = if (modality == "mRNA") generated_data$m else generated_data$g)
  cols <- if (modality == "mRNA") c("dodgerblue4", "deepskyblue1") else c("red", "coral")
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = counts)) +
    ggplot2::geom_histogram(ggplot2::aes(y=..count.., fill = p), binwidth = 1, alpha = 0.7, position = "identity", color = "black") +
    ggplot2::geom_density(alpha = 0.6, adjust = 2) + ggplot2::xlab("UMIs in cell") +
    cowplot::theme_half_open(font_size = 11) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::theme(legend.position = c(0.6, 0.8), legend.title = ggplot2::element_blank()) + ggplot2::labs(title = modality) + ggplot2::scale_fill_manual(values = cols)
}




#' Run EM algo with fast init
#'
#' Runs the EM algorithm using the fast initialization strategy.
#' The assumption is that n is big and pi is small. It is also assumed
#' (for now) that the covariate matrix is non-null.
#'
#' NOTE: Maybe instead take precomputations, i.e. distillation offsets, as arguments.
#'
#'
#' @param dat a data frame containing the columns "m" and "g"
#' @param g_fam family object describing mRNA counts
#' @param m_fam family object describing gRNA counts
#' @param covariate_matrix a data frame storing the covariates
#' @param m_offset the vector of offsets for mRNA model
#' @param g_offset the vector of offsets for gRNA model
#' @param alpha confidence level of CIs
#' @param n_em_rep number of times to repeat EM algorithm on "reduced" data
#'
#' @return fitted model
#' @export
#' @examples
#' m_fam <- g_fam <- augment_family_object(poisson())
#' n <- 200000
#' lib_size <- rpois(n = n, lambda = 10000)
#' m_offset <- g_offset <- log(lib_size)
#' pi <- 0.005
#' m_intercept <- log(0.05)
#' m_perturbation <- log(0.75)
#' g_intercept <- log(0.025)
#' g_perturbation <- log(1.25)
#' covariate_matrix <- data.frame(batch = rbinom(n = n, size = 1, prob = 0.5))
#' m_covariate_coefs <- log(0.9)
#' g_covariate_coefs <- log(1.1)
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept,
#' m_perturbation = m_perturbation, g_fam = g_fam, g_intercept = g_intercept,
#' g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = covariate_matrix, m_covariate_coefs = m_covariate_coefs,
#' g_covariate_coefs = g_covariate_coefs, m_offset = m_offset, g_offset = g_offset)[[1]]
#' m <- dat$m
#' g <- dat$g
#' fit <- run_em_algo_fast_init(m, g, m_fam, g_fam, covariate_matrix, m_offset, g_offset)
run_em_algo_fast_init <- function(m, g, m_fam, g_fam, covariate_matrix, m_offset, g_offset, n_em_rep = 15, pi_guess_range = c(1e-5, 0.02), m_perturbation_guess_range = log(c(0.1, 1.5)), g_perturbation_guess_range = log(c(0.5, 10)), alpha = 0.95) {
  # run the mRNA and gRNA precomputations
  fit_m_precomp <- stats::glm(formula = m ~ ., family = m_fam, data = covariate_matrix, offset = m_offset)
  fit_m_precomp_coefs <- coef(fit_m_precomp)
  fit_g_precomp <- stats::glm(formula = g ~ ., family = g_fam, data = covariate_matrix, offset = g_offset)
  fit_g_precomp_coefs <- coef(fit_g_precomp)

  # obtain the fitted values
  fitted_vals_m_precomp <- as.numeric(stats::fitted.values(fit_m_precomp))
  fitted_vals_g_precomp <- as.numeric(stats::fitted.values(fit_g_precomp))

  # obtain random starting points for reduced GLM-EIV model
  set.seed(4)
  pi_guess <- runif(n = n_em_rep, min = pi_guess_range[1], max = pi_guess_range[2])
  m_perturbation_guess <- runif(n = n_em_rep, min = m_perturbation_guess_range[1], max = m_perturbation_guess_range[2])
  g_perturbation_guess <- runif(n = n_em_rep, min = g_perturbation_guess_range[1], max = g_perturbation_guess_range[2])

  # fit the reduced GLM-EIV model n_em_rep times
  reduced_fits <- lapply(seq(1L, n_em_rep), function(i) {
    run_univariate_poisson_em_algo(m = m, g = g, exp_m_offset = fitted_vals_m_precomp, exp_g_offset = fitted_vals_g_precomp,
                                   m_pert_guess = m_perturbation_guess[i], g_pert_guess = g_perturbation_guess[i], pi_guess = pi_guess[i])
  })

  # among fits that converged, select the one with greatest log-likelihood
  converged <- sapply(reduced_fits, function(fit) fit$converged)
  reduced_fit <- select_best_em_run(reduced_fits[converged])

  # obtain membership probabilities to initialize EM algo
  technical_factors <- colnames(covariate_matrix)
  initial_Ti1s <- run_e_step_pilot(m = m,
                                   g = g,
                                   m_fam = m_fam,
                                   g_fam = g_fam,
                                   pi_guess = reduced_fit$pi,
                                   m_intercept_guess = fit_m_precomp_coefs[["(Intercept)"]],
                                   m_perturbation_guess = reduced_fit$m_perturbation,
                                   m_covariate_coefs_guess = fit_m_precomp_coefs[[technical_factors]],
                                   g_intercept_guess = fit_g_precomp_coefs[["(Intercept)"]],
                                   g_perturbation_guess = reduced_fit$g_perturbation,
                                   g_covariate_coefs_guess = fit_g_precomp_coefs[[technical_factors]],
                                   covariate_matrix = covariate_matrix,
                                   m_offset = m_offset,
                                   g_offset = g_offset)

  # run em algo with initial weights
  fit_em <- run_em_algo_given_weights(m = m, g = g, m_fam = m_fam, g_fam = g_fam, covariate_matrix = covariate_matrix,
                                      initial_Ti1s = initial_Ti1s$Ti1s, m_offset = m_offset, g_offset = g_offset, prev_log_lik = initial_Ti1s$log_lik)
  out <- list(run_inference_on_em_fit(fit_em, alpha = alpha), fit_em)
  return(out)
}

#' Run EM algo given initialization weights
#'
#' Runs GLM-EIV algo with M-step first. Use `run_em_algo_given_pilot` to GLM-EIV algo with E-step first.
#'
#' @param m observed vector m
#' @param g observed vector g
#' @param m_fam augmented family of m
#' @param g_fam augmented family of g
#' @param covariate_matrix the matrix of covariates; set to NULL if there are no covariates.
#' @param initial_Ti1s the initial vector of membership probabilities
#' @param m_offset offsets for GLM for M
#' @param g_offset offsets for GLM for G
#' @param ep_tol (optional) EM convergence threshold
#' @param max_it  (optional) maximum number of EM iterations
#'
#' @return a list containing the following: fit object of M GLM, fit object of G GLM,
#' fit for pi, number of iterations,full log-likelihood of final model
#' @export
#' @name run_em_algo_given_init
#' @examples
#' m_fam <- g_fam <- augment_family_object(poisson())
#' n <- 5000
#' lib_size <- rpois(n = n, lambda = 5000)
#' m_offset <- g_offset <- log(lib_size)
#' pi <- 0.1
#' m_intercept <- log(0.05)
#' m_perturbation <- log(0.8)
#' g_intercept <- log(0.025)
#' g_perturbation <- log(1.2)
#' covariate_matrix <- data.frame(batch = rbinom(n = n, size = 1, prob = 0.5))
#' m_covariate_coefs <- log(0.9)
#' g_covariate_coefs <- log(1.1)
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept,
#' m_perturbation = m_perturbation, g_fam = g_fam, g_intercept = g_intercept,
#' g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = covariate_matrix, m_covariate_coefs = m_covariate_coefs,
#' g_covariate_coefs = g_covariate_coefs, m_offset = m_offset, g_offset = g_offset)[[1]]
#' pi_guess <- 0.05
#' m_intercept_guess <- log(0.07)
#' m_perturbation_guess <- log(0.7)
#' g_intercept_guess <- log(0.02)
#' g_perturbation_guess <- log(1.4)
#' m_covariate_coefs_guess <- log(0.8)
#' g_covariate_coefs_guess <- log(1.2)
#' m <- dat$m
#' g <- dat$g
#' # obtain initial membership probabilities (i.e., run E step) using pilot estimates
#' initial_Ti1s <- run_e_step_pilot(m, g, m_fam, g_fam, pi_guess,
#' m_intercept_guess, m_perturbation_guess, m_covariate_coefs_guess,
#' g_intercept_guess, g_perturbation_guess, g_covariate_coefs_guess,
#' covariate_matrix, m_offset, g_offset)
#' # run em algo
#' fit <- run_em_algo_given_weights(m, g, m_fam, g_fam, covariate_matrix,
#' initial_Ti1s, m_offset, g_offset)
run_em_algo_given_weights <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, prev_log_lik = -Inf, ep_tol = 1e-4, max_it = 75) {
  # augment family objects, if necessary
  if (is.null(m_fam$augmented)) m_fam <- augment_family_object(m_fam)
  if (is.null(g_fam$augmented)) g_fam <- augment_family_object(g_fam)

  # verify column names ok
  check_col_names(covariate_matrix)

  # define some basic quantities
  n <- length(m)
  iteration <- 1L
  converged <- FALSE
  curr_Ti1s <- initial_Ti1s
  prev_log_lik <- prev_log_lik
  log_liks <- numeric()

  # define augmented responses, offsets, and covariate matrix
  augmented_inputs <- augment_inputs(covariate_matrix, m, g, m_offset, g_offset, n)

  # iterate through E and M steps until convergence
  while (!converged) {
    # m step
    m_step <- run_m_step(curr_Ti1s,
                         augmented_inputs$m_augmented, m_fam, augmented_inputs$m_offset_augmented,
                         augmented_inputs$g_augmented, g_fam, augmented_inputs$g_offset_augmented,
                         augmented_inputs$Xtilde_augmented, n)
    curr_log_lik <- m_step$curr_log_lik
    log_liks <- c(log_liks, curr_log_lik)
    curr_tol <- abs(curr_log_lik - prev_log_lik)/min(abs(curr_log_lik), abs(prev_log_lik))
    if (curr_tol < ep_tol) {
      # convergence acheived
      converged <- TRUE
    } else {
      # e step
      curr_Ti1s <- run_e_step(m_step, m, m_fam, g, g_fam, n)
      prev_log_lik <- curr_log_lik
      iteration <- iteration + 1L
      # check iteration limit
      if (iteration >= max_it) {
        break()
      }
    }
  }
  out <- list(fit_m = m_step$fit_m, fit_g = m_step$fit_g, fit_pi = m_step$fit_pi, n_iterations = iteration, log_liks = log_liks, log_lik = curr_log_lik, converged = converged, n = n, posterior_perturbation_probs = m_step$posterior_perturbation_probs)
  return(out)
}


##################
# helper functions
##################
augment_inputs <- function(covariate_matrix, m, g, m_offset, g_offset, n) {
  if (is.null(covariate_matrix)) {
    Xtilde_augmented <- data.frame(perturbation = c(rep(0, n), rep(1, n)))
  } else {
    Xtilde_0 <- dplyr::mutate(covariate_matrix, perturbation = 0)
    Xtilde_1 <- dplyr::mutate(covariate_matrix, perturbation = 1)
    Xtilde_augmented <- rbind(Xtilde_0, Xtilde_1) %>% dplyr::select(perturbation, everything())
  }
  m_augmented <- c(m, m)
  g_augmented <- c(g, g)
  m_offset_augmented <- if (!is.null(m_offset)) c(m_offset, m_offset) else NULL
  g_offset_augmented <- if (!is.null(g_offset)) c(g_offset, g_offset) else NULL
  out <- list(Xtilde_augmented = Xtilde_augmented,
              m_augmented = m_augmented,
              g_augmented = g_augmented,
              m_offset_augmented = m_offset_augmented,
              g_offset_augmented = g_offset_augmented)
  return(out)
}


run_m_step <- function(curr_Ti1s, m_augmented, m_fam, m_offset_augmented, g_augmented, g_fam, g_offset_augmented, Xtilde_augmented, n) {
  # curr_Ti1s[curr_Ti1s < 1e-4] <- 0 # induce weight sparsity
  weights <- c(1 - curr_Ti1s, curr_Ti1s)
  s_curr_Ti1s <- sum(curr_Ti1s)
  fit_pi <- s_curr_Ti1s/n
  if (fit_pi >= 0.5) { # subtract by 1 to ensure label consistency
    s_curr_Ti1s <- n - s_curr_Ti1s
    fit_pi <- s_curr_Ti1s/n
  }
  pi_log_lik <- log(1 - fit_pi) * (n - s_curr_Ti1s) + log(fit_pi) * s_curr_Ti1s

  # fit models for m and g
  m_form <- stats::formula("m_augmented ~ .")
  fit_m <- stats::glm(formula = m_form, data = Xtilde_augmented, family = m_fam,
                      weights = weights, offset = m_offset_augmented)

  g_form <-  stats::formula("g_augmented ~ .")
  fit_g <- stats::glm(formula = g_form, data = Xtilde_augmented, family = g_fam,
                      weights = weights, offset = g_offset_augmented)

  # compute the log-likelihoods
  m_log_lik <- m_fam$get_log_lik(fit_m)
  g_log_lik <- g_fam$get_log_lik(fit_g)

  curr_log_lik <- m_log_lik + g_log_lik + pi_log_lik

  # return list of fitted models, as well as current log-likelihood
  out <- list(fit_m = fit_m, fit_g = fit_g, fit_pi = fit_pi, curr_log_lik = curr_log_lik, posterior_perturbation_probs = curr_Ti1s)
  return(out)
}


update_membership_probs <- function(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi) {
  quotient <- log(1 - fit_pi) - log(fit_pi) + m_fam$d_log_py(m, m_mus_pert0, m_mus_pert1) + g_fam$d_log_py(g, g_mus_pert0, g_mus_pert1)
  out <- 1/(exp(quotient) + 1)
  return(out)
}


# run_e_step <- function(m_step, m, m_fam, g, g_fam, n) {
#  compute_conditional_means <- function(fit, fam, n) {
#    mus <- fam$linkinv(as.numeric(fit$linear.predictors))
#    mus_pert0 <- mus[seq(1, n)]
#    mus_pert1 <- mus[seq(n + 1, 2 * n)]
#    out <- list(mus_pert0 = mus_pert0, mus_pert1 = mus_pert1)
#  }
#  # compute conditional means
#  m_mus <- compute_conditional_means(m_step$fit_m, m_fam, n)
#  g_mus <- compute_conditional_means(m_step$fit_g, g_fam, n)
#  # define all relevant variables
#  fit_pi <- m_step$fit_pi
#  m_mus_pert0 <- m_mus$mus_pert0; m_mus_pert1 <- m_mus$mus_pert1
#  g_mus_pert0 <- g_mus$mus_pert0; g_mus_pert1 <- g_mus$mus_pert1
#  # compute membership probabilities
#  Ti1s <- update_membership_probs(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi)
#  return(Ti1s)
# }


run_e_step_pilot <- function(m, g, m_fam, g_fam, pi_guess, m_intercept_guess, m_perturbation_guess, m_covariate_coefs_guess, g_intercept_guess, g_perturbation_guess, g_covariate_coefs_guess, covariate_matrix, m_offset, g_offset) {
  # compute the conditional means
  m_conditional_means <- compute_theoretical_conditional_means(intercept = m_intercept_guess,
                                                               perturbation_coef = m_perturbation_guess,
                                                               fam = m_fam,
                                                               covariate_matrix = covariate_matrix,
                                                               covariate_coefs = m_covariate_coefs_guess,
                                                               offset = m_offset)
  g_conditional_means <- compute_theoretical_conditional_means(intercept = g_intercept_guess,
                                                               perturbation_coef = g_perturbation_guess,
                                                               fam = g_fam,
                                                               covariate_matrix = covariate_matrix,
                                                               covariate_coefs = g_covariate_coefs_guess,
                                                               offset = g_offset)
  # assign to variables for convenience
  m_mus_pert0 <- m_conditional_means$mu0; m_mus_pert1 <- m_conditional_means$mu1
  g_mus_pert0 <- g_conditional_means$mu0; g_mus_pert1 <- g_conditional_means$mu1
  # compute membership probabilities
  Ti1s <- update_membership_probs(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, pi_guess)
  # compute the log-likelihood
  log_lik <- compute_weighted_log_lik(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, pi_guess, Ti1s)
  return(list(Ti1s = Ti1s, log_lik = log_lik))
}


compute_weighted_log_lik <- function(m_fam, g_fam, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1, fit_pi, Ti1s) {
  # weighted pi log-likelihood; assumes that the Ti1s are correct, i.e., sum(Ti1s)/n < 0.5.
  s_curr_Ti1s <- sum(Ti1s)
  pi_log_lik <- log(1 - fit_pi) * (n - s_curr_Ti1s) + log(fit_pi) * s_curr_Ti1s
  # mRNA log likelihood
  mRNA_log_lik <- m_fam$weighted_log_lik(y = m, mu_0 = m_mus_pert0, mu_1 = m_mus_pert1, Ti1s = Ti1s)
  # gRNA log likelihood
  gRNA_log_lik <- g_fam$weighted_log_lik(y = g, mu_0 = g_mus_pert0, mu_1 = g_mus_pert1, Ti1s = Ti1s)
  # log likelihood
  log_lik <- pi_log_lik + mRNA_log_lik + gRNA_log_lik
  return(log_lik)
}
