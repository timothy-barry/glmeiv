#' Get DFs for mixture plots
#'
#' Returns DFs to plot the mixture distributions.
#'
#' @param sim_spec a simulatr specifier object
#' @param xgrid grid over which to compute the density
#'
#' @return a list containing DFs to plot and optimal thresholds for both mRNA and gRNA
#' @export
get_theoretical_densities_and_thresholds <- function(sim_spec, xgrid) {
  # extract m_pert, m_intercept, g_pert, g_intercept, and pi.
  params <- c("m_intercept", "m_perturbation", "g_intercept", "g_perturbation", "pi")
  param_list <- vector(mode = "list", length = length(params))
  # extract the family objects (assuming these are fixed parameters)
  for (i in seq(1, length(params))) {
    param_list[[i]] <- simulatr::get_param_from_simulatr_spec(simulatr_spec = sim_spec, row_idx = NULL, param = params[i])
  }
  names(param_list) <- params
  # put all the parameters into a data frame
  param_df <- as.data.frame(param_list)
  # extract m_fam, g_fam
  m_fam <- simulatr::get_param_from_simulatr_spec(simulatr_spec = sim_spec, row_idx = NULL, param = "m_fam")
  g_fam <- simulatr::get_param_from_simulatr_spec(simulatr_spec = sim_spec, row_idx = NULL, param = "g_fam")
  # obtain Bayes optimal threshold for each param setting
  n_param_settings <- nrow(param_df)
  # compute plotting data frames
  get_data_frames_to_plot <- function(param_df, fam, modality) {
    fam <- augment_family_object(fam)
    perturbation_name <- paste0(modality, "_perturbation")
    intercept_name <- paste0(modality, "_intercept")
    # compute theoretical conditional means
    conditonal_means <- apply(X = param_df, MARGIN = 1, FUN = function(r) {
      compute_theoretical_conditional_means(intercept = r[[intercept_name]],
                                            perturbation_coef = r[[perturbation_name]],
                                            fam = fam) %>% unlist()
    }) %>% t()
    # obtain the list of density data frames
    plotting_dfs <- lapply(seq(1, nrow(param_df)), function(i) {
      mu0 <- conditonal_means[[i, "mu0"]]
      mu1 <- conditonal_means[[i, "mu1"]]
      pi <- param_df[[i, "pi"]]
      f0 <- fam$density(mu = mu0, xgrid = xgrid) * (1 - pi)
      f1 <- fam$density(mu = mu1, xgrid = xgrid) * pi
      f <- f0 + f1
      rbind(tibble::tibble(f = f0, x = xgrid, component = "Perturbation 0"),
            tibble::tibble(f = f1, x = xgrid, component = "Perturbation 1"),
            tibble::tibble(f = f, x = xgrid, component = "Marginal")) %>%
        dplyr::mutate(component = factor(component, levels = c("Marginal", "Perturbation 0", "Perturbation 1"))) %>%
        dplyr::arrange(component)
    })
    optimal_thresh <- sapply(seq(1, nrow(param_df)), function(i) {
      fam$bayes_classifier(mu0 = conditonal_means[[i, "mu0"]],
                           mu1 =  conditonal_means[[i, "mu1"]],
                           pi = param_df[[i, "pi"]])
      })
    list(plotting_dfs = plotting_dfs, optimal_thresh = optimal_thresh)
  }
  m_out <- get_data_frames_to_plot(param_df, m_fam, "m")
  g_out <- get_data_frames_to_plot(param_df, g_fam, "g")
  out <- list(m_dfs = m_out$plotting_dfs,
       m_thresholds = m_out$optimal_thresh,
       g_dfs = g_out$plotting_dfs,
       g_threshold = g_out$optimal_thresh)
  return(out)
}


#' Plot all arms
#'
#' @param summarized_results summarized results (as outputted by summarize_results) data frame in tibble form.
#' @param parameter parameter to plot
#' @param metric metric to plot
#' @param ylim (optional) ylim
#' @param plot_discont_points (default false) plot vertical lines showing the points where the threshold changes discontinuously? If TRUE, g_thresh must be present as a column in the summarized_results df.
#' @param arm_info (optional) a list containing information about the arms; entries sohuld be "arm_names," "varying_param," and "all_params;" if empty, will be guessed from column names.
#' @param theoretical_values (optional) theoretical values for each arm to plot as horizontal lines
#' @param ylab (optional) y-axis label
#'
#' @return a cowplot containing the plotted arms
#' @export
plot_all_arms <- function(summarized_results, parameter, metric, ylim = NULL, plot_discont_points = FALSE, arm_info = NULL, theoretical_values = NULL, ylab = NULL) {
  if (is.null(arm_info)) {
    arm_info <- list()
    arm_info$arm_names <- grep(pattern = "^arm", x = colnames(summarized_results), value = TRUE)
    arm_info$varying_param <- gsub(pattern = "^arm_", replacement = "", x = arm_info$arm_names)
    arm_info$all_params <- arm_info$varying_param
  }

  summarized_results_sub <- dplyr::filter(summarized_results, parameter == !!parameter, metric == !!metric)
  y_int <- switch(metric, "bias" = 0, "coverage" = 0.95, "count" = NULL, "mse" = 0, "se" = 0)
  arms <- arm_info$arm_names
  ps <- lapply(seq(1, length(arms)), function(i) {
    arm <- arms[i]
    varying_param <- arm_info$varying_param[i]
    fixed_params <- arm_info$all_params[arm_info$all_params != varying_param]
    to_plot <- dplyr::filter(summarized_results_sub, !!as.symbol(arm))
    title <- sapply(fixed_params, function(fixed_param) paste0(fixed_param, " = ",
                                                               to_plot[[fixed_param]][1])) %>% paste0(., collapse = ", ")
    if (plot_discont_points) {
      to_plot_thresh <- to_plot %>% dplyr::filter(method == "thresholding") %>% dplyr::arrange(arm)
      g_thresh_floor <- to_plot_thresh %>% dplyr::pull(g_thresh) %>% floor()
      diff_thresh <- c(diff(g_thresh_floor), 0)
      xintercepts <- to_plot_thresh[which(diff_thresh != 0), varying_param] %>% dplyr::pull()
    }
    theoretical_values_passed <- !is.null(theoretical_values)
    if (theoretical_values_passed) {
      theoretical_value_sub <- dplyr::filter(theoretical_values, arm == !!arm)
    }

    p <- ggplot2::ggplot(to_plot, ggplot2::aes(x = !!as.symbol(varying_param), y = value, col = method)) +
      ggplot2::geom_hline(yintercept = y_int, lwd = 0.6) +
      (if (plot_discont_points) ggplot2::geom_vline(xintercept = xintercepts, lwd = 0.3, col = "lightslategray")) +
      ggplot2::geom_point() + (if (!theoretical_values_passed) ggplot2::geom_line() else NULL) + (ggplot2::ylab( if (is.null(ylab)) metric else ylab[i] )) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower_mc_ci, ymax = upper_mc_ci), width = 0) +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5)) +
      ggplot2::ggtitle(title) + (if (is.null(ylim)) NULL else ggplot2::ylim(ylim)) +
      ggplot2::scale_color_manual(values = c("dodgerblue3", "firebrick1"), breaks = c("em", "thresholding")) +
      (if (theoretical_values_passed) ggplot2::geom_line(data = theoretical_value_sub, mapping = ggplot2::aes(x = x, y = value), linetype = "dashed") else NULL)


    l <- cowplot::get_legend(p + ggplot2::theme(legend.position = "bottom"))
    p_out <- p + ggplot2::theme(legend.position = "none")
    return(list(plot = p_out, legend = l, n_method = to_plot$method %>% unique() %>% length()))
  })
  n_ps <- length(ps)
  legend_idx <- which.max(sapply(ps, function(i) i$n_method))
  vert_plot <- cowplot::plot_grid(plotlist = lapply(ps, function(i) i$plot),
                                  align = "v", axis = "l", nrow = n_ps, labels = letters[1:n_ps])
  out <- cowplot::plot_grid(vert_plot, ps[[legend_idx]]$legend, ncol = 1, rel_heights = c(1, .1))
  return(out)
}


#' Plot mixture
#'
#' Plots a mixture distribution.
#'
#' @param density_df a density df as outputted by "get_theoretical_densities_and_thresholds."
#' @param thresh a threshold to plot as a vertical line
#' @param x_max maximal x-value
#' @param x_min minimal x-value
#' @param points plot points? (if FALSE, only lines)
#' @param xlab x-axis label
#'
#' @return a ggplot object
#' @export
plot_mixture <- function(density_df, thresh = NULL, x_max = Inf, x_min = -Inf, points = TRUE, xlab = NULL) {
  if (!is.null(x_max)) density_df <- dplyr::filter(density_df, x <= x_max, x >= x_min)
  p <- ggplot2::ggplot(data = density_df, mapping = ggplot2::aes(x = x, y = f, col = component)) + (if (!is.null(thresh)) ggplot2::geom_vline(xintercept = thresh, lwd = 0.3) else NULL) + ggplot2::geom_line() + (if (points) ggplot2::geom_point() else NULL) + ggplot2::theme_bw() + ggplot2::ylab("Density") + (if (!is.null(xlab)) ggplot2::xlab(xlab) else NULL)
  return(p)
}


#' Plot posterior membership probabilities
#'
#' Plots all sets of posterior membership probabilities for a given grid_row_id of a sim_res data frame.
#'
#' @param sim_res raw result data frame outputted by simulatr
#' @param grid_row_id an grid row id
#' @param valid_ids character vector of valid IDs
#' @param n_approx_01 location on y-axis at which to draw horizontal line
#'
#' @return a list of length 2: (i) a list of histogram plots, and (ii) a data frame storing summary metrics for each run
#' @export
plot_posterior_membership_probabilities <- function(sim_res, grid_row_id, valid_ids, n_approx_01 = 50) {
  sim_res_w_mem_probs <- sim_res %>% dplyr::filter(method == "em", grid_row_id == !!grid_row_id) %>% dplyr::group_by(id) %>% dplyr::filter("membership_prob" %in% target)
  all_runs <- sim_res_w_mem_probs %>% dplyr::group_map(.f = function(tbl, key) {
    p <- ggplot2::ggplot(data = dplyr::filter(tbl, target == "membership_prob"),
                         mapping = ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(binwidth = 0.1, col = "black", fill = "gray") + ggplot2::theme_bw() +
      ggplot2::xlab("Membership probability") + ggplot2::ylab("Count") + ggplot2::geom_hline(yintercept = n_approx_01, col = "darkred") +
      ggplot2::geom_vline(xintercept = 0.15, col = "darkblue") + ggplot2::geom_vline(xintercept = 0.85, col = "darkblue")
    id <- as.character(key %>% dplyr::pull())
    return(list(p = p, id = id))
  })
  metrics_df <- sim_res_w_mem_probs %>% dplyr::summarize(spread = value[target == "membership_probability_spread"],
                                                         n_approx_0 = value[target == "n_approx_0"],
                                                         n_approx_1 = value[target == "n_approx_1"],
                                                         m_perturbation = value[parameter == "m_perturbation" & target == "estimate"],
                                                         g_perturbation = value[parameter == "g_perturbation" & target == "estimate"],
                                                         pi = value[parameter == "pi" & target == "estimate"],
                                                         valid = id[1] %in% valid_ids) %>% dplyr::mutate(id = as.character(id))
  plots <- lapply(all_runs, function(run) run$p)
  names(plots) <- sapply(all_runs, function(run) run$id)
  return(list(plots = plots, metrics_df = metrics_df))
}


#' Plot EM classifications
#'
#' Creates a histogram containing the EM estimates for m_perturbation colored by type.
#'
#' There are four categories: the cartesian product of (confident vs. unconfident) and (pluasible vs. implausible).
#' An estimate is "confident" if the posterior membership probabilities are sufficiently spread out (see function
#' obtain valid IDs); and estimate is "plausible" if the estimates themselves are in biologically realistic ranges.
#' Valid IDs are those that are both confident and plausible.
#'
#' @param em_classifications dataframe outputted by `obtain_valid_ids`
#' @param sim_spec a simulatr specifier object (to supply the ground truth)
#' @param grid_row_id grid_row_id giving row to plot
#' @param categories classification categories to include in plot
#' @param parameter parameter to make histogram for
#'
#' @return a histogram (in ggplot format)
#' @export
plot_em_classifications <- function(em_classifications, sim_spec, grid_row_id, categories = c("confident-plausible", "confident-implausible", "unconfident-plausible", "unconfident-implausible"), parameter = "m_perturbation") {
  em_classifications_to_plot <- dplyr::filter(em_classifications, grid_row_id == !!grid_row_id, classification %in% categories)
  gt <- simulatr::get_param_from_simulatr_spec(sim_spec, grid_row_id, param = parameter)
  cols <- c("confident-plausible" = "firebrick4", "confident-implausible" = "maroon4", "unconfident-plausible" = "dodgerblue4", "unconfident-implausible" = "lightseagreen")
  p <- ggplot2::ggplot(data = em_classifications_to_plot, mapping = ggplot2::aes(x = !!as.symbol(parameter), col = classification)) + ggplot2::geom_histogram(bins = 25, alpha = 0.7, position = "identity", fill = "gray") + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = gt) + ggplot2::scale_color_manual(values = cols) + ggplot2::geom_hline(yintercept = 0)
  return(p)
}


classify_estimates_em <- function(sim_res, spread_thresh = 0.1, approx_0_thresh = 60, approx_1_thresh = 60, g_pert_lower = -0.25, g_pert_upper = Inf, m_pert_lower = -Inf, m_pert_upper = 0.25, pi_lower = 0, pi_upper = 0.3) {
  out <- dplyr::filter(sim_res, method == "em") %>% dplyr::group_by(id) %>%
    dplyr::summarize(
      confident_output = (value[target == "converged"] == 1
                          & value[target == "membership_probability_spread"] > spread_thresh
                          & value[target == "n_approx_0"] >= approx_0_thresh
                          & value[target == "n_approx_1"] >= approx_1_thresh),
      g_perturbation = value[parameter == "g_perturbation" & target == "estimate"],
      m_perturbation = value[parameter == "m_perturbation" & target == "estimate"],
      pi = value[parameter == "pi" & target == "estimate"],
      plausible_estimates = (g_perturbation >= g_pert_lower &
                               g_perturbation <= g_pert_upper &
                               m_perturbation >= m_pert_lower &
                               m_perturbation <= m_pert_upper &
                               pi >= pi_lower & pi <= pi_upper),
      grid_row_id = grid_row_id[1]) %>%
    dplyr::mutate(classification = paste0(ifelse(confident_output, "confident", "unconfident"), "-",
                                          ifelse(plausible_estimates, "plausible", "implausible")) %>% factor(),
                  valid = confident_output & plausible_estimates)
  return(out)
}


#' Obtain valid IDs
#'
#' Obtains valid IDs given the output of a simulatr experiment.
#'
#' @param sim_res sim_res object
#' @param spread_thresh spread threshold
#' @param approx_0_thresh number cells approximately 0 threshold
#' @param approx_1_thresh number cells approximately 1 threshold
#' @param g_pert_lower minimum value for g_pert
#' @param g_pert_upper maximum value for g_pert
#' @param m_pert_lower minimum value for m_pert
#' @param m_pert_upper maximum value for m_pert
#' @param pi_lower minimum value for pi
#' @param pi_upper maximum value for pi
#'
#' @return a list of length two; (i) a data frame giving the classification of each EM run, and (ii) a character vector of valid IDs for both EM and thresholding
#' @export
obtain_valid_ids <- function(sim_res, spread_thresh = 0.1, approx_0_thresh = 50, approx_1_thresh = 50, g_pert_lower = -0.3, g_pert_upper = Inf, m_pert_lower = -Inf, m_pert_upper = 0.3, pi_lower = 0, pi_upper = 0.5) {
  em_df <- classify_estimates_em(sim_res, spread_thresh, approx_0_thresh, approx_1_thresh, g_pert_lower, g_pert_upper, m_pert_lower, m_pert_upper, pi_lower, pi_upper)
  valid_thresh_ids <- sim_res %>%
    dplyr::filter(method == "thresholding") %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(valid_idx = (value[target == "fit_attempted"] == 1)) %>%
    dplyr::filter(valid_idx) %>% dplyr::pull(id) %>% as.character()
  valid_em_ids <- dplyr::filter(em_df, valid) %>% dplyr::pull(id) %>% as.character()
  valid_ids <- c(valid_em_ids, valid_thresh_ids)
  return(list(em_classifications = em_df, valid_ids = valid_ids))
}


plot_profile_likelihood <- function(param_to_set, val_to_set, fit, m_offsets, g_offsets, covariate_matrix, m_fam, g_fam) {
  # m coefs
  m_coefs <- fit$fit_m$coefficients
  m_perturbation <- as.list(m_coefs)[["perturbation"]]
  m_intercept <- as.list(m_coefs)[["(Intercept)"]]
  m_covariate_coefs <- if (length(m_coefs) >= 3) as.list(m_coefs)[[seq(3, length(m_coefs))]] else NULL

  # g coefs
  g_coefs <- fit$fit_g$coefficients
  g_perturbation <- as.list(g_coefs)[["perturbation"]]
  g_intercept <- as.list(g_coefs)[["(Intercept)"]]
  g_covariate_coefs <- if (length(g_coefs) >= 3) as.list(g_coefs)[[seq(3, length(g_coefs))]] else NULL

  # pi
  pi <- fit$fit_pi

  # update value of given parameter
  assign(x = param_to_set, value = val_to_set)

  # pi, m, g
  m <- as.numeric(fit$fit_m$y)
  g <- as.numeric(fit$fit_g$y)

  # compute conditional means
  m_conditional_means <- compute_theoretical_conditional_means(intercept = m_intercept,
                                                               perturbation_coef = m_perturbation,
                                                               fam = m_fam,
                                                               covariate_matrix = covariate_matrix,
                                                               covariate_coefs = m_covariate_coefs,
                                                               offset = m_offsets)
  g_conditional_means <- compute_theoretical_conditional_means(intercept = g_intercept,
                                                               perturbation_coef = g_perturbation,
                                                               fam = g_fam,
                                                               covariate_matrix = covariate_matrix,
                                                               covariate_coefs = g_covariate_coefs,
                                                               offset = g_offsets)
  # run e step to get fitted probabilities and log likelihood
  e_step <- run_e_step(m_fam = m_fam, g_fam = g_fam, m = m, g = g,
             m_mus_pert0 = m_conditional_means$mu0, m_mus_pert1 = m_conditional_means$mu1,
             g_mus_pert0 = g_conditional_means$mu0, g_mus_pert1 = g_conditional_means$mu1,
             fit_pi = pi)
  return(e_step$log_lik)
}
