#' @title Plot selection biased joint density of sample size and effect size
#' @description Plot selection biased joint density of sample size and effect size
#' @param n Vector of sample sizes
#' @param d Vector of effect sizes
#' @param delta Effect size mean
#' @param exp_delta Expected effect size
#' @param sigma2 Effect size variance
#' @param n_size Negative binomial size parameter
#' @param n_mu Negative binomial mean parameter
#' @param alpha Type I error rate
#' @param beta Type II error rate
#' @param p_pbs weights of probability  selection given p(select | p > alpha)
#' @param p_ssa weights of probability selection given p(select | n < n_req)
#' @return ggplot2 object
#' @export
#' @examples
#' plot_joint_density(
#'   delta = 0.4, exp_delta = 0.4, sigma2 = 0.3, n_size = 4,
#'   n_mu = 120, p_pbs = 0.5, p_ssa = 0.5
#'   )

plot_joint_density <- function(n = seq(10, 250, 1), d = seq(-2, 2, 0.01), delta, exp_delta, sigma2, n_size, n_mu, alpha = 0.05, beta = 0.2, p_pbs, p_ssa) {
  requireNamespace("ggtext", quietly = TRUE)
  requireNamespace("cetcolor", quietly = TRUE)
  requireNamespace("scales", quietly = TRUE)
  n_lims <- range(n)
  d_lims <- range(d)
  nreq <- pwr::pwr.t.test(d = exp_delta, power = 1 - beta)$n
  grid_data <- expand.grid(d = d, n = n)
  grid_data$density <- stats::dnorm(grid_data$d, mean = delta, sd = sqrt(sigma2)) * stats::dnbinom(grid_data$n, size = n_size, mu = n_mu)
  grid_data$p <- p_from_nd(grid_data$n, grid_data$d)
  grid_data$w_pbs <- ifelse(grid_data$p < 0.05, 1, p_pbs)
  grid_data$w_ssa <- ifelse(grid_data$n > nreq, 1, p_ssa)
  grid_data$density_weighted <- with(grid_data, w_pbs*w_ssa*density)

  p <- ggplot2::ggplot(grid_data, ggplot2::aes(n, d)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::geom_raster(ggplot2::aes(fill = density_weighted)) +
    ggplot2::geom_linerange(
      ggplot2::aes(x = ceiling(nreq), ymin = d_lims[1], ymax = d_from_np(nreq, TRUE)),
      color = "white", linewidth = 0.1
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(x = ceiling(nreq), ymin = d_lims[2], ymax = d_from_np(nreq, FALSE)),
      color = "white", linewidth = 0.1
    ) +
    ggplot2::geom_line(
      data = data.frame(n = seq(n_lims[1], n_lims[2], length.out = 1e4)),
      ggplot2::aes(x = n, y = d_from_np(n, TRUE)),
      color = "white",
      linewidth = 0.1,
      linetype = "solid",
      lineend = "square"
    ) +
    ggplot2::geom_line(
      data = data.frame(n = seq(n_lims[1], n_lims[2], length.out = 1e4)),
      ggplot2::aes(x = n, y = d_from_np(n, FALSE)),
      color = "white",
      linetype = "solid",
      lineend = "square",
      linewidth = 0.1
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "white", linewidth = .1, linetype = "longdash") +
    ggplot2::scale_fill_gradientn(
      colors = cetcolor::cet_pal(100, "l8"),
      name = "Dichte",
      breaks = scales::breaks_pretty(n = 4),
      limits = c(0, NA)
    ) +
    ggplot2::scale_x_continuous(
      name = "Stichprobengr&ouml;sse *N*",
      limits = n_lims,
      breaks = pretty(n_lims, n = 6),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      name = "Effektst&auml;rke<br>Cohens *d*",
      limits = d_lims,
      breaks = pretty(d_lims, n = 6),
      expand = c(0, 0)
    ) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggtext::element_markdown(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggtext::element_markdown(margin = ggplot2::margin(r = 10), lineheight = 1.25),
      plot.background = ggplot2::element_blank(),
      legend.position = "left",
      legend.title = ggplot2::element_text(colour = "grey30"),
      legend.text = ggplot2::element_text(colour = "grey30"),
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm"),
      axis.line = ggplot2::element_line(colour = "grey80"),
      axis.text = ggplot2::element_text(colour = "grey50"),
      axis.title = ggplot2::element_text(colour = "grey30")
    )
  return(suppressWarnings(print(p)))
  }

