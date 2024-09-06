#' @export
print.speec_control <- function(x, ...) {
  c_bounds <- Map(function(lower, upper) c(lower, upper), x$bounds$lower, x$bounds$upper)
  format_bounds <- sapply(names(c_bounds), function(var) {
    paste0(var, ": [", toString(c_bounds[[var]]), "]")
  }, USE.NAMES = FALSE)

  not_specified_start <- sapply(x$start, is.null)
  x$start[not_specified_start] <- "automatic, see doc."

  format_start <- sapply(names(x$start), function(var) paste0(var, ": [", toString(x$start[[var]]), "]"))
  if (is.null(x$hyperparameters$vf)) x$hyperparameters$vf <- "not specified"

  cli::cli({
  cli::cli_h1("Control Parameters for SPEEC Publication Bias Estimation")
  cli::cli_h2("KDE Estimation Parameters:")
  cli::cli_li(
    items = c(
      "{.field bw} - Bandwidth method: {.val {x$bw}}",
      "{.field n_grid} - Number of grid points: {.val {x$n_grid}}",
      "{.field pr} - CDF quantiles: {.val {x$pr}}"
      ),
  )
  cli::cli_h2("Simulation parameters:")
  cli::cli_li(
    items = c(
      "{.field k_sim} - Number of simulated primary studies: {.val {x$k_sim}}",
      "{.field alpha} - Type I error rate: {.val {x$alpha}}",
      "{.field beta} - Type II error rate: {.val {x$beta}}",
      "{.field slope_ssp} - Slope for the weighting function of sample size planning: {.val {x$slope_ssp}}"
    )
  )
  cli::cli_h2("Optimization parameters:")
  cli::cli_li(
    items = c(
      "{.field bounds} ({.fn set_boundaries}):  {.val {format_bounds}}",
      "{.field start} ({.fn set_start}) - Starting values: {format_start}"
    )
  )
  cli::cli_h2("Simulated Annealing Optimization Parameters")
  cli::cli_text("Via {.fn set_hyperparameters}:")
  cli::cli_li(
    items = c(
      "{.field vf} - Variation function: {.val {x$hyperparameters$vf}}",
      "{.field rf} - Random factor: {.val {x$hyperparameters$rf}}",
      "{.field dyn_rf} - Dynamic random factor: {.val {x$hyperparameters$dyn_rf}}",
      "{.field t0} - Initial temperature: {.val {x$hyperparameters$t0}}",
      "{.field nlimits} - Maximum number of iterations: {.val {x$hyperparameters$nlimit}}",
      "{.field r} - Temperature reduction (outer loop): {.val {x$hyperparameters$r}}",
      "{.field k} - Constant for the Metropolis function: {.val {x$hyperparameters$k}}",
      "{.field t_min} - Minimum temperature (stops outer loop): {.val {x$hyperparameters$t_min}}",
      "{.field maxgood} - Maximum number of loss function improvements (stops inner loop): {.val {x$hyperparameters$maxgood}}",
      "{.field stopac} - Maximum number of repetitions with low loss improvement (stops inner loop): {.val {x$hyperparameters$stopac}}",
      "{.field ac_acc} - Accuracy of the stopac break criterion: {.val {x$hyperparameters$ac_acc}}"
    )
  )
  })
  invisible(x)
}

#' @export
plot.speec_optim <- function(x, type = c("loss", "rf", "par")) {
    trace_data <- x$trace
    if(length(trace_data) == 0) {
        cli::cli_abort("Trace argument not set to {.val TRUE}. Can´t plot the results.")
    } else {
        type <- rlang::arg_match(type)
        if(type == "loss") {
            p <- ggplot2::ggplot(
                data = trace_data,
                mapping = ggplot2::aes(x = n_outer, y = loss)
            ) +
            ggplot2::labs(x = "Outer loop iterations", y = "Loss")
        } else if(type == "rf") {
            trace_data <- trace_data[, grep("(^rf)|(n_outer)", names(trace_data))]
            trace_data_long <- stats::reshape(
                data = trace_data,
                varying = list(c("rf_phi_n", "rf_mu_n", "rf_mu_d", "rf_sigma2_d", "rf_w_pbs")),
                v.names = "value",
                timevar = "parameter",
                times = c("rf_phi_n", "rf_mu_n", "rf_mu_d", "rf_sigma2_d", "rf_w_pbs"),
                direction = "long"
                )
            p <- ggplot2::ggplot(
                data = trace_data_long,
                mapping = ggplot2::aes(x = n_outer, y = value)
            ) +
            ggplot2::facet_wrap(~parameter, scales = "free_y") +
            ggplot2::labs(x = "Outer loop iterations", "Random factors")
        } else if(type == "par") {
            trace_data <- trace_data[, grep("(_n$)|(_d$)|(n_outer)|(w_pbs)", names(trace_data))]
            trace_data <- trace_data[, !grepl("^rf_", names(trace_data))]
            trace_data_long <- stats::reshape(
                data = trace_data,
                varying = list(c("phi_n", "mu_n", "mu_d", "sigma2_d", "w_pbs")),
                v.names = "value",
                timevar = "parameter",
                times = c("phi_n", "mu_n", "mu_d", "sigma2_d", "w_pbs"),
                direction = "long"
                )
            p <- ggplot2::ggplot(
                data = trace_data_long,
                mapping = ggplot2::aes(x = n_outer, y = value)
            ) +
            ggplot2::facet_wrap(~parameter, scales = "free_y") +
            ggplot2::labs(x = "Outer loop iterations", "Parameter")
        }
        p <- p +
            ggplot2::geom_line() +
            ggplot2::theme_bw()
        return(p)
    }
}
