
simetabias_control <- function(
    bw = c("silverman", "scott", "sheather-jones", "ucv", "bcv"),
    n_grid = c(2**7+1, 2**7+1), quantile = c(0.005, 0.995), k = 1e4,
    alpha = .05, beta = .2, discr = 4, d_delta_hat_bounds = c(0, 3), w_pbs_bounds = c(0, 1), d_delta_hat_tol = .1,
    alpha_bounds = .001, start_w_pbs = .5, trace = TRUE, only_pbs = FALSE
) {
  if(length(n_grid) != 2) cli::cli_abort("n_grid must be a vector of length 2.")
  if(any(quantile < 0) || any(quantile > 1)) cli::cli_abort("quantile must be between 0 and 1.")
  if(length(quantile) != 2) cli::cli_abort("quantile must be a vector of length 2.")
  if(quantile[1] >= quantile[2]) cli::cli_abort("the first element must be smaller than the second element.")
  if(length(d_delta_hat_bounds) != 2) cli::cli_abort("d_delta_hat_bounds must be a vector of length 2.")
  if(length(w_pbs_bounds) != 2) cli::cli_abort("w_pbs_bounds must be a vector of length 2.")
  if(w_pbs_bounds[1] >= w_pbs_bounds[2]) cli::cli_abort("the first element must be smaller than the second element.")
  if(w_pbs_bounds[1] < 0 || w_pbs_bounds[2] > 1) cli::cli_abort("w_pbs_bounds must be between 0 and 1.")
  if(d_delta_hat_bounds[1] >= d_delta_hat_bounds[2]) cli::cli_abort("the first element must be smaller than the second element.")
  bw <- rlang::arg_match(bw)
  if (only_pbs) discr <- 0
  out <- named_list(bw, n_grid, quantile, k, alpha, beta, discr, d_delta_hat_bounds,
             w_pbs_bounds, d_delta_hat_tol, alpha_bounds, start_w_pbs, trace,
             only_pbs)
  class(out) <- "simetabias_control"
  return(out)
}

sa_control <- function(vf = NULL, rf = 1, dyn_rf = TRUE, t0 = 1e3, nlimit = 100,
                       r = 0.6, k = 1, t_min = 0.1, maxgood = 100, stopac = 30, ac_acc = 1e-4) {
  named_list(vf, rf, dyn_rf, t0, nlimit, r, k, t_min, maxgood, stopac, ac_acc)
}

# get_boundaries <- function(emp_data, d_delta_hat_bounds = c(0, 3), w_pbs_bounds = c(0, 1), alpha_bounds = .001) {
#   # for n
#   n <- emp_data[, "n"]
#   fit <- MASS::glm.nb(n ~ 1, link = log)
#   ci_nb_mu <- exp(coef(fit) + sqrt(diag(vcov(fit))) * qnorm(c(alpha_bounds/2, 1 - alpha_bounds/2)))
#   if(ci_nb_mu[1] < 0) ci_nb_mu[1] <- 0
#   ci_nb_phi <- fit$theta + fit$SE.theta * qnorm(c(alpha_bounds/2, 1 - alpha_bounds/2))
#   if(ci_nb_phi[1] < 0) ci_nb_phi[1] <- 0
#
#   # for d
#   d <- emp_data[, "d"]
#   k <- length(d)
#   ci_norm_sd <- sqrt((k - 1) * var(d) / qchisq(c(1 - alpha_bounds/2, alpha_bounds/2), df = k - 1))
#   ci_norm_mean <- mean(d) + sd(d)/sqrt(k) * qnorm(c(alpha_bounds/2, 1 - alpha_bounds/2))
#
#   L <- list(ci_nb_phi, ci_nb_mu, ci_norm_mean, ci_norm_sd, d_delta_hat_bounds, w_pbs_bounds)
#   out <- list(
#     lower = sapply(L, "[[", 1),
#     upper = sapply(L, "[[", 2)
#   )
#   return(out)
# }

get_boundaries <- function(
    n_phi_bounds = c(1, 50), n_mu_bounds = c(50, 300),
    d_mu_bounds = c(0, 1), d_sigma_bounds = c(0.1, 1),
    d_delta_hat_bounds = c(0, 3), w_pbs_bounds = c(0, 1),
    only_pbs = FALSE) {
  if (only_pbs) {
    L <- list(n_phi_bounds, n_mu_bounds, d_mu_bounds, d_sigma_bounds, w_pbs_bounds)
  } else {
    L <- list(n_phi_bounds, n_mu_bounds, d_mu_bounds, d_sigma_bounds, d_delta_hat_bounds, w_pbs_bounds)
  }
  out <- list(
    lower = sapply(L, "[[", 1),
    upper = sapply(L, "[[", 2)
  )
  return(out)
}

get_start <- function(emp_data, d_delta_hat_tol = .1, start_w_pbs = 0.5, only_pbs = FALSE) {
  start_nb <- estimate_nb(emp_data[, "n"])
  start_norm <- estimate_norm(emp_data[, "d"])
  start_d_delta_hat <- mean(emp_data[, "d"])
  if (start_d_delta_hat == 0) start_d_delta_hat <- d_delta_hat_tol
  if (only_pbs) {
    out <- c(start_nb, start_norm, w_pbs = start_w_pbs)
  } else {
    out <- c(start_nb, start_norm, start_d_delta_hat = start_d_delta_hat, w_pbs = start_w_pbs)
  }
  return(out)
}

optimize_sa <- function(emp_data, simetabias_control = simetabias_control(), sa_control = sa_control()) {
  bw <- simetabias_control$bw
  n_grid <- simetabias_control$n_grid
  q <- simetabias_control$quantile
  k <- simetabias_control$k
  alpha <- simetabias_control$alpha
  beta <- simetabias_control$beta
  discr <- simetabias_control$discr
  d_delta_hat_bounds <- simetabias_control$d_delta_hat_bounds
  d_delta_hat_tol <- simetabias_control$d_delta_hat_tol
  w_pbs_bounds <- simetabias_control$w_pbs_bounds
  start_w_pbs <- simetabias_control$start_w_pbs
  trace <- simetabias_control$trace
  ac_acc <- sa_control$ac_acc
  alpha_bounds <- simetabias_control$alpha_bounds
  only_pbs <- simetabias_control$only_pbs

  lims <- find_kde_limits(emp_data, quantile = q)
  fhat_empirical <- bivariate_kde(emp_data, bw = bw, n_grid = n_grid, lims = lims)

  if (only_pbs) {
    f <- function(params) {
      theo_data <- simulate_meta(
        k = k, n_size = params[1], n_mu = params[2], d_mu = params[3], d_sigma = params[4],
        d_delta_hat = 1, w_pbs = params[5], discr = discr, beta = beta, alpha = alpha
      )
      fhat_theoretical <- bivariate_kde(theo_data, bw = bw, n_grid = n_grid, lims = lims)
      kl_div(fhat_empirical, fhat_theoretical)
    }
  } else {
    f <- function(params) {
      theo_data <- simulate_meta(
        k = k, n_size = params[1], n_mu = params[2], d_mu = params[3], d_sigma = params[4],
        d_delta_hat = params[5], w_pbs = params[6], discr = discr, beta = beta, alpha = alpha
      )
      fhat_theoretical <- bivariate_kde(theo_data, bw = bw, n_grid = n_grid, lims = lims)
      kl_div(fhat_empirical, fhat_theoretical)
    }
  }

  start_params <- get_start(emp_data, d_delta_hat_tol = d_delta_hat_tol, start_w_pbs = start_w_pbs, only_pbs = only_pbs)
  bounds <- get_boundaries(d_delta_hat_bounds = d_delta_hat_bounds, w_pbs_bounds = w_pbs_bounds, only_pbs = only_pbs)

  sa_control$ac_acc <- ac_acc*f(start_params)
  opt <- optimization::optim_sa(
    fun = f, start = start_params,
    lower = bounds$lower, upper = bounds$upper,
    maximization = FALSE,
    trace = trace,
    control = sa_control
  )
  return(opt)
}




