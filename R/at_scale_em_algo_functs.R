#' Run GLM-EIV given precomputations
#'
#' Runs a full GLM-EIV given mRNA and gRNA precomputations, as outputted by `run_glmeiv_precomputation` function.
#'
#' @param m_precomp mRNA precomputation list, as outputted by `run_glmeiv_precomputation`
#' @param g_precomp gRNA precomputation list, as outputted by `run_glmeiv_precomputation`
#' @param fam_str string describing family; currently "poisson" or "Negative Binomial" supported
#' @param n_em_rep number of EM reps to run on reduced GLM-EIV model
#' @param pi_guess_range range of values over which to sample pi
#' @param m_perturbation_guess_range range of values over which to sample m_perturbation
#' @param g_perturbation_guess_range range of values over which to sample g_perturbation
#' @inheritParams run_full_glmeiv_given_weights
#'
#' @return a fitted GLM-EIV model
#' @export
run_glmeiv_given_precomputations <- function(m, g, m_precomp, g_precomp, covariate_matrix, m_offset, g_offset, n_em_rep, pi_guess_range, m_perturbation_guess_range, g_perturbation_guess_range) {
  set.seed(4)
  m_fam <- m_precomp$fam
  g_fam <- g_precomp$fam

  # get random starting guesses
  guesses <- lapply(list(pi = pi_guess_range,
                         m_perturbation = m_perturbation_guess_range,
                         g_perturbation = g_perturbation_guess_range), function(r) {
                           stats::runif(n = n_em_rep, min = r[1], max = r[2])
                         })

  # fit the reduced GLM-EIV model over the random starting vectors
  exp_m_offset <- m_precomp$fitted_values; exp_g_offset <- g_precomp$fitted_values
  reduced_fits <- lapply(seq(1L, n_em_rep), function(i) {
    run_reduced_em_algo(m = m, g = g, exp_m_offset = exp_m_offset, exp_g_offset = exp_g_offset,
                        m_pert_guess = guesses$m_perturbation[i],
                        g_pert_guess = guesses$g_perturbation[i],
                        pi_guess = guesses$pi[i], m_fam = m_fam, g_fam = g_fam)
    })

  # obtain the best run according to log-likelihood
  best_run <- select_best_em_run(reduced_fits)

  # Finally, run full GLM-EIV using pilot estimates
  fit <- run_full_glmeiv_given_pilot_params(m = m, g = g, m_fam = m_fam, g_fam = g_fam,
                                            pi_guess = best_run$pi, m_intercept_guess = m_precomp$fitted_intercept,
                                            m_perturbation_guess = best_run$m_perturbation, m_covariate_coefs_guess = m_precomp$covariate_coefs,
                                            g_intercept_guess = g_precomp$fitted_intercept, g_perturbation_guess = best_run$g_perturbation,
                                            g_covariate_coefs_guess = g_precomp$covariate_coefs, covariate_matrix = covariate_matrix,
                                            m_offset = m_offset, g_offset = g_offset)
  return(fit)
}


#' Run GLM-EIV precomputation
#'
#' Runs a precomputation on a gene or gRNA.
#'
#' @param y a vector of gene or gRNA expressions
#' @param covariate_matrix the matrix of technical factors
#' @param offset a vector of offsets, generally log library size
#' @param fam_str a string specifying the family, currently either "Negative Binomial" or "poisson"
#'
#' @return a list containing (i) estimated theta (in NB case only), (ii) fitted values, (iii) fitted intercept term, (iv) covariate_coef terms
#' @export
run_glmeiv_precomputation <- function(y, covariate_matrix, offset, fam) {
  out <- list()
  fam_str <- fam$fam_str
  if (fam_str == "Negative Binomial" && is.na(fam$theta)) {
    form <- stats::as.formula(if (is.null(offset)) "y ~ ." else "y ~ . + offset(offset)")
    fit_precomp <- MASS::glm.nb(formula = form, data = dplyr::mutate(data.frame(y = y), covariate_matrix))
    fam <- augment_family_object(MASS::negative.binomial(fit_precomp$theta))
  } else if (fam_str %in% c("poisson", "Negative Binomial")) {
    fit_precomp <- stats::glm(formula = y ~ ., family = fam, data = dplyr::mutate(data.frame(y = y), covariate_matrix), offset = offset)
  } else stop("Family string not recognized.")
  coefs <- as.list(stats::coef(fit_precomp))
  fitted_values <- as.numeric(stats::fitted.values(fit_precomp))
  fitted_intercept <- coefs[["(Intercept)"]]
  covarite_coefs <- unlist(coefs[names(coefs) != "(Intercept)"])
  out$fitted_values <- fitted_values
  out$fitted_intercept <- fitted_intercept
  out$covariate_coefs <- covarite_coefs
  out$fam <- fam
  return(out)
}


#' Run GLM-EIV at scale simulatr
#'
#' @param dat data frame containing columns "m" and "g" for mRNA counts and gRNA counts
#' @param alpha confidence level
#' @param n_em_rep number of times to repeat EM algorithm in reduced model
#' @param save_membership_probs_mult save posterior membership probabilities at this multiple
#' @param fam_str string (currently either "poisson" or "Negative Binomial") giving family object
#' @return fitted GLM-EIV model
#' @export
#' @inheritParams run_full_glmeiv_given_weights
#' @inheritParams run_glmeiv_given_precomputations
#'
#' @examples
#' m_fam <- g_fam <- augment_family_object(poisson())
#' n <- 200000
#' lib_size <- rpois(n = n, lambda = 5000)
#' m_offset <- g_offset <- log(lib_size)
#' pi <- 0.005
#' m_intercept <- log(0.01)
#' m_perturbation <- log(0.5)
#' g_intercept <- log(0.005)
#' g_perturbation <- log(2.5)
#' covariate_matrix <- data.frame(batch = rbinom(n = n, size = 1, prob = 0.5))
#' m_covariate_coefs <- log(0.9)
#' g_covariate_coefs <- log(1.1)
#' dat <- generate_full_data(m_fam = m_fam, m_intercept = m_intercept,
#' m_perturbation = m_perturbation, g_fam = g_fam, g_intercept = g_intercept,
#' g_perturbation = g_perturbation, pi = pi, n = n, B = 2,
#' covariate_matrix = covariate_matrix, m_covariate_coefs = m_covariate_coefs,
#' g_covariate_coefs = g_covariate_coefs, m_offset = m_offset, g_offset = g_offset)[[1]]
#' m <- dat$m; g <- dat$g; p <- dat$p
#' # ability to recover ground truth given p
#' fit <- run_glmeiv_at_scale_simulatr(dat, m_fam, g_fam, covariate_matrix, m_offset, g_offset)
run_glmeiv_at_scale_simulatr <- function(dat, m_fam, g_fam, covariate_matrix, m_offset, g_offset, alpha = 0.95, n_em_rep = 15, save_membership_probs_mult = 250, pi_guess_range = c(1e-5, 0.03), m_perturbation_guess_range = log(c(0.1, 1.5)), g_perturbation_guess_range = log(c(0.5, 10))) {
  m <- dat$m
  g <- dat$g
  # if precomp already complete...
  if ("m_precomp" %in% names(attributes(dat)) && "g_precomp" %in% names(attributes(dat))) {
    m_precomp <- attr(dat, "m_precomp")
    g_precomp <- attr(dat, "g_precomp")
  } else { # run the precomputations
    m_precomp <- run_glmeiv_precomputation(y = m, covariate_matrix = covariate_matrix, offset = m_offset, fam = m_fam)
    g_precomp <- run_glmeiv_precomputation(y = g, covariate_matrix = covariate_matrix, offset = g_offset, fam = g_fam)
  }
  # run glmeiv given precomputations, timing it.
  time <- system.time({
    fit <- run_glmeiv_given_precomputations(m = m, g = g, m_precomp = m_precomp, g_precomp = g_precomp, covariate_matrix = covariate_matrix, m_offset = m_offset, g_offset = g_offset, n_em_rep = n_em_rep, pi_guess_range = pi_guess_range, m_perturbation_guess_range = m_perturbation_guess_range, g_perturbation_guess_range = g_perturbation_guess_range)
    s <- run_inference_on_em_fit(fit, alpha)
    })[["elapsed"]]
  # do post-processing (by a call to a function), then return result.
  out <- wrangle_glmeiv_result(s, time, fit)
  return(out)
}


#' Run GLM-EIV (random inititalization)
#'
#' Runs GLM-EIV using random initializations for the parameters. Currently, the function assumes that there is (at most) a single covariate term.
#'
#' @param alpha CIs returned at level (1-alpha)%
#' @param n_em_rep number of EM replicates
#' @param save_membership_probs_mult save membership probabilities at this multiple
#' @param m_intercept_guess_range range over which to sample m_intercept
#' @param g_intercept_guess_range range over which to sample g_intercept
#' @param m_covariate_coefs_guess_range range over which to sample m_covariate_coefs
#' @param g_covariate_coefs_guess_range range over which to sample g_covariate_coefs
#' @param dat a data frame containing columns m, g
#'
#' @inheritParams run_full_glmeiv_given_weights
#' @inheritParams run_glmeiv_given_precomputations
#'
#' @return a fitted GLM-EIV object
#' @export
run_glmeiv_random_init_simulatr <- function(dat, m_fam, g_fam, covariate_matrix, m_offset, g_offset, alpha = 0.95, n_em_rep = 15, save_membership_probs_mult = 250, pi_guess_range = c(1e-5, 0.03), m_perturbation_guess_range = log(c(0.1, 1.5)), g_perturbation_guess_range = log(c(0.5, 10)), m_intercept_guess_range = log(c(1e-4, 1e-1)), g_intercept_guess_range = log(c(1e-4, 1e-1)), m_covariate_coefs_guess_range = log(c(0.5, 1.5)), g_covariate_coefs_guess_range = log(c(0.5, 1.5))) {
  # get random starting guesses for the parameters; first, five core model parameters
  m <- dat$m
  g <- dat$g
  # if covariate matrix is NULL, set g_covariate_coefs_guess_range, m_covariate_coefs_guess_range to NULL
  if (is.null(covariate_matrix)) {
    g_covariate_coefs_guess_range <- m_covariate_coefs_guess_range <- NULL
  }
  # set the family objects; if available as attribute of dat, use that; else, use the m_fam and g_fam passed as arguments
  if ("m_precomp" %in% names(attributes(dat)) && "g_precomp" %in% names(attributes(dat))) {
    m_precomp <- attr(dat, "m_precomp"); m_fam <- m_precomp$fam
    g_precomp <- attr(dat, "g_precomp"); g_fam <- g_precomp$fam
  }
  best_fit <- NULL
  best_log_lik <- -Inf
  # run method, timing the result
  time <- system.time({
  # generate guesses
  guesses <- lapply(list(pi = pi_guess_range,
                         m_perturbation = m_perturbation_guess_range,
                         g_perturbation = g_perturbation_guess_range,
                         m_intercept = m_intercept_guess_range,
                         g_intercept = g_intercept_guess_range,
                         m_covariate_coefs = m_covariate_coefs_guess_range,
                         g_covariate_coefs = g_covariate_coefs_guess_range), function(r) {
                           if (is.null(r)) NULL else stats::runif(n = n_em_rep, min = r[1], max = r[2])
                         })
    # fit models for each starting guess
    for (i in seq(1L, n_em_rep)) {
      print(paste0("Rep ", i))
      fit <- run_full_glmeiv_given_pilot_params(m = m, g = g, m_fam = m_fam, g_fam = g_fam,
                                                pi_guess = guesses$pi[i],
                                                m_intercept_guess = guesses$m_intercept[i],
                                                m_perturbation_guess = guesses$m_perturbation[i],
                                                m_covariate_coefs_guess = guesses$m_covariate_coefs[i],
                                                g_intercept_guess = guesses$g_intercept[i],
                                                g_perturbation_guess = guesses$g_perturbation[i],
                                                g_covariate_coefs_guess = guesses$g_covariate_coefs[i],
                                                covariate_matrix = covariate_matrix,
                                                m_offset = m_offset, g_offset = g_offset)
      if (fit$log_lik > best_log_lik) {
        best_fit <- fit
        best_log_lik <- fit$log_lik
      }
    }
  s <- run_inference_on_em_fit(best_fit)
  })[["elapsed"]]
  # process the result
  out <- process_glmeiv_results_simulatr(best_fit, s, dat, save_membership_probs_mult, time)
  return(out)
}


#' Wrangle GLM-EIV result
#'
#' converts output of GLM-EIV into useful form.
#'
#' @param s output of `run_inference_on_em_fit`
#' @param time execution time
#' @param em_fit a fitted GLM-EIV object
#'
#' @return data frame with results in tidy form
#' @export
wrangle_glmeiv_result <- function(s, time, em_fit) {
  s_trans <- dplyr::filter(s, parameter != "pi") %>%
    dplyr::select(parameter, estimate, p_value, confint_lower, confint_upper) %>%
    dplyr::mutate(estimate = exp(estimate),
                  confint_lower = exp(confint_lower),
                  confint_upper = exp(confint_upper)) %>%
    dplyr::add_row(dplyr::filter(s, parameter == "pi") %>%
                     dplyr::select(parameter, estimate, p_value, confint_lower, confint_upper))

  membership_prob_spread <- compute_mean_distance_from_half(em_fit$posterior_perturbation_probs)
  n_approx_1 <- sum(em_fit$posterior_perturbation_probs > 0.85)
  n_approx_0 <- sum(em_fit$posterior_perturbation_probs < 0.15)
  # transform output result
  meta_df <- tibble::tibble(parameter = "meta",
                            target = c("converged", "membership_probability_spread",
                                       "n_approx_0", "n_approx_1", "time"),
                            value = c(em_fit$converged, membership_prob_spread, n_approx_0, n_approx_1, time))
  out <- rbind(tidyr::pivot_longer(s_trans, cols = -parameter, names_to = "target"), meta_df)
  return(out)
}


#' Select best EM run
#'
#' Returns the best run, defined as the run with m_perturbation in a reasonable range and
#' greatest log-likelihood.
#'
#' @param em_runs a list of outputs of run_em_algo_given_init
#' @param m_perturbation_range range of acceptable values for m_perturbation
#'
#' @return the best run of the list
#' @export
select_best_em_run <- function(em_runs, m_perturbation_range = c(log(0.1), log(2))) {
  # restrict runs to those with m_perturbation in acceptable range
  m_perturbations <- sapply(em_runs, function(run) run$m_perturbation)
  valid_m_pert <- (m_perturbations >= m_perturbation_range[1]) & (m_perturbations <= m_perturbation_range[2])
  # if there exists at least one run with valid m_pert, subset
  if (any(valid_m_pert)) em_runs <- em_runs[valid_m_pert]
  # return the em run with greatest log-likelihood
  log_liks <- sapply(em_runs, function(run) run$log_lik)
  return(em_runs[[which.max(log_liks)]])
}
