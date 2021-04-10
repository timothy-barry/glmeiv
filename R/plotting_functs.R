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


plot_em_runs <- function(em_runs) {
  to_plot <- lapply(seq(1, length(em_runs)), function(i) {
    curr_log_liks <- em_runs[[i]]$log_liks
    data.frame(log_lik = curr_log_liks, iteration = seq(1, length(curr_log_liks))) %>%
      dplyr::mutate(run = i)
    }) %>% do.call(what = rbind, args = .) %>% dplyr::mutate(run = factor(run))
  ggplot2::ggplot(to_plot %>% dplyr::filter(iteration >= 3), ggplot2::aes(x = iteration, y = log_lik, col = run)) +
    ggplot2::geom_line() + ggplot2::theme_bw()
}
