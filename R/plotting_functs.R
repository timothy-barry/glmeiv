#' Get DFs for mixture plots
#'
#' Returns DFs to plot the mixture distributions.
#'
#' @param sim_spec a simulatr specifier object
#' @param x_grid grid over which to compute the density
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
#'
#' @return
plot_all_arms <- function(summarized_results, parameter, metric, ylim = NULL, plot_discont_points = FALSE, arm_info = NULL, theoretical_values = NULL) {
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
      ggplot2::geom_point() + (if (!theoretical_values_passed) ggplot2::geom_line() else NULL) + ggplot2::ylab(metric) +
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
#'
#' @return a ggplot object
#' @export
plot_mixture <- function(density_df, thresh = NULL, x_max = NULL, x_min = NULL, points = TRUE) {
  if (!is.null(x_max)) density_df <- dplyr::filter(density_df, x <= x_max, x >= x_min)
  p <- ggplot2::ggplot(data = density_df, mapping = ggplot2::aes(x = x, y = f, col = component)) + (if (!is.null(thresh)) ggplot2::geom_vline(xintercept = thresh, lwd = 0.3) else NULL) + ggplot2::geom_line() + (if (points) ggplot2::geom_point() else NULL) + ggplot2::theme_bw() + ggplot2::ylab("Density")
  return(p)
}
