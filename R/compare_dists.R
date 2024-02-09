compare_control <- function(emp_dist = NULL, kde_bounds = c(0.005, 0.995), n_rep = 1, k = 1e3, n_grid = 100,
                            start_nb = c(10, 100), start_norm = c(0, 0.5),
                            discr = 5, alpha = 0.05, beta = 0.2, parallel = FALSE, workers = NULL) {

  sa_list <- named_list(
    emp_dist, kde_bounds,n_rep,k,n_grid, start_nb, start_norm,
    discr, alpha, beta, parallel, workers
    )
  if (is.null(emp_dist)) cli::cli_abort("No empirical distribution provided.")
  return(sa_list)
}

compare_dists <- function(params, control = compare_control()) {
  if (control$n_rep > 1) {
    if (control$parallel) {
      requireNamespace("future", quietly = TRUE)
      requireNamespace("future.apply", quietly = TRUE)
      if (is.null(control$workers)) {
        workers <- future::availableCores() - 1
        future::plan(future::multisession, workers = workers)
      } else {
        if (!is.numeric(control$workers)) cli::cli_abort("Number of workers must be numeric.")
        if (control$workers > future::availableCores()) cli::cli_abort("Number of workers is larger than available cores.")
        if (control$workers <= 0) cli::cli_abort("Number of workers must be positive.")
        future::plan(future::multisession, workers = control$workers)
      }
      dist_vec <- future.apply::future_replicate(
        n = control$n_rep,
        expr = {
          theo_dist <- simulate_meta(
            k = control$k, k_size = params[1], k_mu = params[2], delta = params[3], sigma2 = params[4],
            d_delta_hat = params[5], w_pbs = params[6],
            discr = control$discr, beta = control$beta, alpha = control$alpha
          )
          compute_loss(control$emp_dist, theo_dist, control$n_grid, control$kde_bounds, control$start_nb, control$start_norm)
        },
        simplify = TRUE,
      )
    } else {
      dist_vec <- replicate(
        n = control$n_rep,
        expr = {
          theo_dist <- simulate_meta(
            k = control$k, k_size = params[1], k_mu = params[2], delta = params[3], sigma2 = params[4],
            d_delta_hat = params[5], w_pbs = params[6],
            discr = control$discr, beta = control$beta, alpha = control$alpha
          )
          compute_loss(control$emp_dist, theo_dist, control$n_grid, control$kde_bounds, control$start_nb, control$start_norm)
        },
        simplify = TRUE
      )
    }
    tvd <- mean(dist_vec)
  } else if (control$n_rep == 1) {
    theo_dist <- simulate_meta(
      k = control$k, k_size = params[1], k_mu = params[2], delta = params[3], sigma2 = params[4],
      d_delta_hat = params[5], w_pbs = params[6],
      discr = control$discr, beta = control$beta, alpha = control$alpha
    )
    tvd <- compute_loss(control$emp_dist, theo_dist, control$n_grid, control$kde_bounds, control$start_nb, control$start_norm)
  } else {
    cli::cli_abort("Number of replications must be positive.")
  }
  return(tvd)
}
