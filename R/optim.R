#' @title Control of Important SPEEC Parameters
#' @description
#' Contruct control structures for the SPEEC publication bias estimation framework
#' @param bw bandwidth selection method. One of "silverman" (\link[stats]{bw.nrd0}), "scott" (\link[stats]{bw.nrd}),
#' "sheather-jones" (\link[stats]{bw.SJ}), "ucv" (\link[stats]{bw.ucv}) or "bcv" (\link[stats]{bw.bcv}).
#' @param n_grid integer vector of length 2, specifying the number evenly spaced grid points along each axis for density estimation.
#' @param pr numeric vector of length 2, specifying the lower and upper quantile for the CDF
#' @param k_sim integer vector of length 1, specifying the number of simulated primary studies for each iteration of the optimization algorithm
#' @param alpha fixed type I error rate
#' @param beta fixed type II error rate
#' @param slope_ssp integer vector of length 1, specifying the discrimination parameter for the weighting function of sample size planning
#' @param bounds a list generated with \code{\link{set_boundaries}} specifying the lower and upper bounds for the optimization algorithm
#' @param only_pbs Logical scalar. If TRUE, the optimization algorithm will only consider the publication bias selection parameter. Default is FALSE.
#' @param trace Logical scalar. If TRUE, interim results are stored. Necessary for the plot function. Default is TRUE.
#' @param hyperparameters a list generated with \code{\link{set_hyperparameters}} specifying the hyperparameters for the simulated annealing optimization algorithm
#' @export
#'
speec_control <- function(bw = c("silverman", "scott", "sheather-jones", "ucv", "bcv"),
                          n_grid = c(2**7 + 1, 2**7 + 1), pr = c(0.005, 0.995), k_sim = 1e4,
                          bounds = set_boundaries(), alpha = .05, beta = .2, slope_ssp = 4, only_pbs = TRUE,
                          hyperparameters = set_hyperparameters()
                          ) {
  if (length(n_grid) != 2) cli::cli_abort("n_grid must be a vector of length 2.")
  if (any(pr < 0) || any(pr > 1)) cli::cli_abort("pr must be between 0 and 1.")
  if (length(pr) != 2) cli::cli_abort("{.arg pr} must be a vector of length 2.")
  if (pr[1] >= pr[2]) cli::cli_abort("the first element in {.arg pr} must be smaller than the second element.")
  if (!rlang::is_scalar_integerish(k_sim)) cli::cli_abort("{.arg i} must be a scalar integer.")
  if (k_sim < 5000) cli::cli_alert_warning("{.arg k_sim} is less than a 5000. This is not recommended!")
  bw <- rlang::arg_match(bw)
  if (only_pbs) {
    cli::cli_alert_info("Running the optimization algorithm with {.code only_pbs = TRUE}.")
    cli::cli_alert_info("The delta_hat parameter will be ignored and {.arg slope_ssa} will be set to 0")
    slope_ssp <- 0
    bounds$lower <- bounds$lower[which(names(bounds$lower) != "delta_hat")]
    bounds$upper <- bounds$upper[which(names(bounds$upper) != "delta_hat")]
  }
  out <- named_list(bw, n_grid, pr, k_sim, alpha, beta, slope_ssp, only_pbs, bounds, hyperparameters)
  class(out) <- "speec_control"
  return(out)
}

#' @title Control of Important Hyperparameters for the Differential Evolution Optimization
#' @description
#' @VTR 
#' Contruct control structures for the differential evolution optimization algorithm
#' @export
#'
set_hyperparameters <- function(
  VTR = -Inf, strategy = 1, NP = 150, itermax = 1000,
  CR = 0.1, F = 0.9, trace = TRUE, initialpop = NULL, storepopfrom = 1, 
  storepopfreq = 1, p = 0.2, c = 0, bs = TRUE, reltol = sqrt(.Machine$double.eps),
  steptol = itermax
  ) {
    RcppDE::DEoptim.control(
      VTR = VTR, strategy = strategy, bs = bs, NP = NP,
      itermax = itermax, CR = CR, F = F, trace = trace,
      initialpop = initialpop, storepopfrom = storepopfrom,
      storepopfreq = storepopfreq, p = p, c = c, reltol = reltol,
      steptol = steptol
    )
}

#' Boundary Constraints for Simulated Annealing Optimization
#' @description
#' Construct a list of lower and upper boundaries for the optimization algorithm.
#' @param phi_n integer vector of length 2, specifying the lower and upper boundaries for the dispersion parameter of the negative binomial distribution
#' @param mu_n integer vector of length 2, specifying the lower and upper boundaries for the mean of the negative binomial distribution
#' @param mu_d numeric vector of length 2, specifying the lower and upper boundaries for the mean of the normal distribution
#' @param sigma2_d numeric vector of length 2, specifying the lower and upper boundaries for the standard deviation of the normal distribution
#' @param delta_hat numeric vector of length 2, specifying the lower and upper boundaries for the expected effect size for sample size planning
#' @param w_pbs numeric vector of length 2, specifying the lower and upper boundaries for the publiation bias weights
#' @return
#' A list with two elements, lower and upper, each containing a named list of the lower and upper boundaries for the optimization algorithm.
#' @export
set_boundaries <- function(
    phi_n = c(1, 50), mu_n = c(50, 500), mu_d = c(-2, 2), sigma2_d = c(0.1, 2),
    delta_hat = c(0.01, 3), w_pbs = c(0, 1)
) {
  check_arg_len(2, phi_n, mu_n, mu_d, sigma2_d, delta_hat, w_pbs)
  L <- named_list(phi_n, mu_n, mu_d, sigma2_d, delta_hat, w_pbs)
  out <- list(
    lower = sapply(L, "[[", 1),
    upper = sapply(L, "[[", 2)
  )
  return(out)
}


#' @title Estimate and Correct Publication Bias in Meta-Analysis under Consideration of Sample Size Planning
#' @param emp_data A dataframe or matrix with two columns: effect size (d) and sample size (n)
#' @param speec_control Control parameters \code{\link{speec_control}} for details.
#' @returns
#' The output is a speec_optim list object with following entries:
#'   \describe{
#'     \item{\code{par}}{
#'       Named parameter values after optimization.
#'     }
#'     \item{\code{function_value}}{
#'       Loss function value after optimization.
#'     }
#'     \item{\code{start}}{
#'       The initial function variables. Set with \code{\link{set_start}}.
#'     }
#'     \item{\code{lower}}{
#'       The lower boundaries of the function variables. Set with \code{\link{set_boundaries}}.
#'     }
#'     \item{\code{upper}}{
#'       The upper boundaries of the function variables. Set with \code{\link{set_boundaries}}.
#'     }
#'     \item{\code{control}}{
#'       Control arguments, see \code{\link{speec_control}} and \code{\link{set_hyperparameters}}.
#'     }
#'     \item{\code{runtime}}{
#'       The runtime of the optimization in seconds.
#'     }
#'   }
#' @export
#' @example man/examples/speec_example.R
speec <- function(emp_data, speec_control = speec_control()) {
  # extract control parameters
  k_sim <- speec_control$k_sim
  bw <- speec_control$bw
  n_grid <- speec_control$n_grid
  alpha <- speec_control$alpha
  beta <- speec_control$beta
  slope_ssp <- speec_control$slope_ssp
  pr <- speec_control$pr
  bounds <- speec_control$bounds
  hyperparameters <- speec_control$hyperparameters
  only_pbs <- speec_control$only_pbs

  # compute empirical kernel density estimate
  lims <- find_kde_limits(emp_data, pr = pr)
  fhat_empirical <- bivariate_kde(emp_data, bw = bw, n_grid = n_grid, lims = lims)

  # optimization function
  if (only_pbs) {
    loss_function <- function(params) {
      theo_data <- simulate_meta(
        k_sim = k_sim, phi_n = params[1], mu_n = params[2], mu_d = params[3], sigma2_d = params[4],
        delta_hat = 1, w_pbs = params[5], slope_ssp = slope_ssp, beta = beta, alpha = alpha, only_pbs = only_pbs
      )
      fhat_theoretical <- bivariate_kde(theo_data, bw = bw, n_grid = n_grid, lims = lims)
      kl_div(fhat_empirical, fhat_theoretical)
    }
  } else {
    loss_function <- function(params) {
      theo_data <- simulate_meta(
        k_sim = k_sim, phi_n = params[1], mu_n = params[2], mu_d = params[3], sigma2_d = params[4],
        delta_hat = params[5], w_pbs = params[6], slope_ssp = slope_ssp, beta = beta, alpha = alpha, only_pbs = only_pbs
      )
      fhat_theoretical <- bivariate_kde(theo_data, bw = bw, n_grid = n_grid, lims = lims)
      kl_div(fhat_empirical, fhat_theoretical)
    }
  }
  t1 <- Sys.time()
  opt <- RcppDE::DEoptim(fn = loss_function, lower = bounds$lower, upper = bounds$upper, control, hyperparameters)
  t2 <- Sys.time()
  opt <- unclass(opt)
  opt$runtime <- difftime(t2, t1, units = "secs")
  class(opt) <- "speec_optim"
  return(opt)
}
