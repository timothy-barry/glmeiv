#' Run EM algo given initialization weights
#'
#' @param m observed vector m
#' @param g observed vector g
#' @param m_fam augmented family of m
#' @param g_fam augmented family of g
#' @param covariate_matrix the matrix of covariates; set to NULL if there are no covariates.
#' @param initial_Ti1s the initial vector of membership probabilities
#' @param m_offset (optional) offsets for GLM for M
#' @param g_offset (optional) offsets for GLM for G
#' @param ep_tol EM convergence threshold
#' @param max_it maximum number of EM iterations
#'
#' @return a list containing the following: fit object of M GLM, fit object of G GLM,
#' fit for pi, number of iterations,full log-likelihood of final model
#' @export
#' @examples
#' n <- 1000
#' m_fam <- augment_family_object(poisson())
#' g_fam <- augment_family_object(poisson())
#' m_offset <- g_offset <- NULL
#' pi <- 0.2
#' covariate_matrix <- data.frame(p_mito = runif(n = n, 0, 10),
#'                                lib_size = rpois(n = n, lambda = 1))
#' m_coef <- c(-1, -2, 1, 0.5)
#' g_coef <- c(-1, 2, 1, 0.5)
#' generated_data <- generate_data_from_model(m_fam, g_fam, m_coef, g_coef, pi, covariate_matrix)
#' m <- generated_data$m
#' g <- generated_data$g
#' initial_Ti1s <- (g >= median(g))
#' em_res <- run_em_algo_given_init(m, g, m_fam, g_fam,
#' covariate_matrix, initial_Ti1s, m_offset, g_offset)
run_em_algo_given_init <- function(m, g, m_fam, g_fam, covariate_matrix, initial_Ti1s, m_offset, g_offset, ep_tol = 0.1, max_it = 50) {
  # verify column names ok
  check_col_names(covariate_matrix)

  # define some basic quantities
  n <- length(m)
  iteration <- 0L
  converged <- FALSE
  curr_Ti1s <- initial_Ti1s
  prev_log_lik <- -Inf

  # define augmented responses, offsets, and covariate matrix
  augmented_inputs <- augment_inputs(covariate_matrix, m, g, m_offset, g_offset)

  # iterate through E and M steps until convergence
  while (!converged) {
    # e step
    e_step <- run_e_step(curr_Ti1s,
                         augmented_inputs$m_augmented, m_fam, augmented_inputs$m_offset_augmented,
                         augmented_inputs$g_augmented, g_fam, augmented_inputs$g_offset_augmented,
                         augmented_inputs$Xtilde_augmented, n)
    curr_log_lik <- e_step$curr_log_lik

    if (abs(curr_log_lik - prev_log_lik) < ep_tol) {
      # convergence acheived
      converged <- TRUE
    } else {
      # m step
      curr_Ti1s <- run_m_step(e_step, m, m_fam, g, g_fam, n)
      prev_log_lik <- curr_log_lik
      iteration <- iteration + 1L
      # check iteration limit
      if (iteration >= max_it) {
        warning("Iteration limit reached; solution not optimal.")
        break()
      }
    }
  }
  out <- list(fit_m = e_step$fit_m, fit_g = e_step$fit_g, fit_pi = e_step$fit_pi, n_iterations = iteration, log_lik = curr_log_lik, converged = converged)
  return(out)
}


##################
# helper functions
##################

augment_inputs <- function(covariate_matrix, m, g, m_offset, g_offset) {
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


run_e_step <- function(curr_Ti1s, m_augmented, m_fam, m_offset_augmented, g_augmented, g_fam, g_offset_augmented, Xtilde_augmented, n) {
  weights <- c(1 - curr_Ti1s, curr_Ti1s)
  # fit augmented models
  fit_m <- stats::glm(formula = m_augmented ~ ., data = Xtilde_augmented, family = m_fam,
                      weights = weights, offset = m_offset_augmented)
  fit_g <- stats::glm(formula = g_augmented ~ ., data = Xtilde_augmented, family = g_fam,
                      weights = weights, offset = g_offset_augmented)
  fit_pi <- sum(curr_Ti1s)/n

  # compute the log-likelihoods
  m_log_lik <- stats::logLik(fit_m)[1]
  g_log_lik <- stats::logLik(fit_g)[1]
  pi_log_lik <- log(1 - fit_pi) * (n - sum(curr_Ti1s)) + log(fit_pi) * sum(curr_Ti1s)
  curr_log_lik <- m_log_lik + g_log_lik + pi_log_lik

  # return list of fitted models, as well as current log-likelihood
  out <- list(fit_m = fit_m, fit_g = fit_g, fit_pi = fit_pi, curr_log_lik = curr_log_lik)
  return(out)
}


run_m_step <- function(e_step, m, m_fam, g, g_fam, n) {
  # compute conditional means
  m_mus <- compute_conditional_means(e_step$fit_m, m_fam, n)
  g_mus <- compute_conditional_means(e_step$fit_g, g_fam, n)
  # define all relevant variables
  fit_pi <- e_step$fit_pi
  m_mus_pert0 <- m_mus$mus_pert0; m_mus_pert1 <- m_mus$mus_pert1
  g_mus_pert0 <- g_mus$mus_pert0; g_mus_pert1 <- g_mus$mus_pert1
  # compute membership probabilities
  update_membership_probs <- update_membership_probs_factory(m_fam, g_fam, fit_pi)
  Ti1s <- mapply(update_membership_probs, m, g, m_mus_pert0, m_mus_pert1, g_mus_pert0, g_mus_pert1)
  return(Ti1s)
}


compute_conditional_means <- function(fit, fam, n) {
  mus <- fam$linkinv(as.numeric(fit$linear.predictors))
  mus_pert0 <- mus[seq(1, n)]
  mus_pert1 <- mus[seq(n + 1, 2 * n)]
  out <- list(mus_pert0 = mus_pert0, mus_pert1 = mus_pert1)
}


update_membership_probs_factory <- function(m_fam, g_fam, fit_pi) {
  m_log_py_given_mu <- m_fam$log_py_given_mu
  g_log_py_given_mu <- g_fam$log_py_given_mu
  f <- function(m_i, g_i, m_mus_i0, m_mus_i1, g_mus_i0, g_mus_i1) {
    alpha_i0 <- m_log_py_given_mu(m_i, m_mus_i0) + g_log_py_given_mu(g_i, g_mus_i0) + 1 - fit_pi
    alpha_i1 <-  m_log_py_given_mu(m_i, m_mus_i1) + g_log_py_given_mu(g_i, g_mus_i1) + fit_pi
    return(1/(1 + exp(alpha_i0 - alpha_i1)))
  }
  return(f)
}
