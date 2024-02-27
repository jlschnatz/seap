#' @export
print.simetabias_control <- function(x, ...) {
  cli::cli_h1("Control Parameters for Simetabias Publication Bias Estimation")
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
  c_bounds <- Map(function(lower, upper) c(lower, upper), x$bounds$lower, x$bounds$upper)
  format_bounds <- sapply(names(c_bounds), function(var) {
    paste0(var, ": [", toString(c_bounds[[var]]), "]")
  }, USE.NAMES = FALSE)

  not_specified_start <- sapply(x$start, is.null)
  x$start[not_specified_start] <- "automatic, see doc."

  format_start <- sapply(names(x$start), function(var) {
    paste0(var, ": [", toString(x$start[[var]]), "]")
  })

  cli::cli_li(
    items = c(
      "{.field bounds} ({.fn set_boundaries}):  {.val {format_bounds}}",
      "{.field start} ({.fn set_start}) - Starting values: {format_start}"
    )
  )
  cli::cli_h2("Simulated Annealing Optimization Parameters")
  if (is.null(x$hyperparameters$vf)) x$hyperparameters$vf <- "not specified"
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
  invisible(x)
}

#' @title Control of Important Simetabias Parameters
#' @description
#' Contruct control structures for the simetabias publication bias estimation framework
#' @param bw bandwidth selection method. One of "silverman" (\link[stats]{bw.nrd0}), "scott" (\link[stats]{bw.nrd}),
#' "sheather-jones" (\link[stats]{bw.SJ}), "ucv" (\link[stats]{bw.ucv}) or "bcv" (\link[stats]{bw.bcv}).
#' @param n_grid integer vector of length 2, specifying the number evenly spaced grid points along each axis for density estimation.
#' @param pr numeric vector of length 2, specifying the lower and upper quantile for the CDF
#' @param k_sim integer vector of length 1, specifying the number of simulated primary studies for each iteration of the optimization algorithm
#' @param alpha fixed type I error rate
#' @param beta fixed type II error rate
#' @param slope_ssp integer vector of length 1, specifying the discrimination parameter for the weighting function of sample size planning
#' @param bounds a list generated with \code{\link{set_boundaries}} specifying the lower and upper bounds for the optimization algorithm
#' @param start a named list generated with \code{\link{set_start}} specifying the starting values for the optimization algorithm
#' @param only_pbs Logical scalar. If TRUE, the optimization algorithm will only consider the publication bias selection parameter. Default is FALSE.
#' @param trace Logical scalar. If TRUE, interim results are stored. Necessary for the plot function. Default is TRUE.
#' @param hyperparameters a list generated with \code{\link{set_hyperparameters}} specifying the hyperparameters for the simulated annealing optimization algorithm
#' @export
#'
simetabias_control <- function(bw = c("silverman", "scott", "sheather-jones", "ucv", "bcv"),
                               n_grid = c(2**7 + 1, 2**7 + 1), pr = c(0.005, 0.995), k_sim = 1e4,
                               bounds = set_boundaries(), start = set_start(),
                               alpha = .05, beta = .2, slope_ssp = 4, only_pbs = FALSE,
                               trace = TRUE, hyperparameters = set_hyperparameters()
                               ) {
  if (length(n_grid) != 2) cli::cli_abort("n_grid must be a vector of length 2.")
  if (any(pr < 0) || any(pr > 1)) cli::cli_abort("pr must be between 0 and 1.")
  if (length(pr) != 2) cli::cli_abort("{.arg pr} must be a vector of length 2.")
  if (pr[1] >= pr[2]) cli::cli_abort("the first element in {.arg pr} must be smaller than the second element.")
  if (!rlang::is_scalar_integerish(k_sim)) cli::cli_abort("{.arg i} must be a scalar integer.")
  if (k_sim < 5000) cli::cli_alert_warning("{.arg i} is less than a 5000. This is not recommended!")
  bw <- rlang::arg_match(bw)
  if (only_pbs) {
    cli::cli_alert_info("Running the optimization algorithm with {.code only_pbs = TRUE}.")
    cli::cli_alert_info("The delta_hat parameter will be ignored and {.arg slope_ssa} will be set to 0")
    slope_ssp <- 0
    bounds$lower <- bounds$lower[which(names(bounds$lower) != "delta_hat")]
    bounds$upper <- bounds$upper[which(names(bounds$upper) != "delta_hat")]
    start <- start[which(names(start) != "delta_hat")]
  }

  out <- named_list(
    bw, n_grid, pr, k_sim, alpha, beta, slope_ssp, trace, only_pbs, bounds, start,
     hyperparameters
  )
  class(out) <- "simetabias_control"
  return(out)
}

#' @title Control of Important Hyperparameters for the Simulated Annealing Optimization
#' @description
#' Contruct control structures for the simulated annealing optimization algorithm
#' @param vf Function that determines the variation of the function variables for the next iteration. The variation function is allowed to depend on the vector of variables of the current iteration, the vector of random factors rf and the temperature of the current iteration. Default is a uniform distributed random number with relative range rf.
#' @param rf Numeric vector. Random factor vector that determines the variation of the random number of vf in relation to the dimension of the function variables for the following iteration. Default is 1. If dyn_rf is enabled, the rf change dynamically over time.
#' @param dyn_rf Logical scalar. rf change dynamically over time to ensure increasing precision with increasing number of iterations. Default is TRUE, see 'details'.
#' @param t0 Numeric scalar. Initial temperature. Default is 1000.
#' @param nlimit Integer scalar. Maximum number of iterations of the inner loop. Default is 100.
#' @param r Numeric. Temperature reduction in the outer loop. Default is 0.6.
#' @param k Numeric. Constant for the Metropolis function. Default is 1.
#' @param t_min Numeric. Temperature where outer loop stops. Default is 0.1.
#' @param maxgood Integer. Break criterion to improve the algorithm performance. Maximum number of loss function improvements in the inner loop. Breaks the inner loop. Default is 100.
#' @param stopac Integer. Break criterion to improve the algorithm performance. Maximum number of repetitions where the loss improvement is lower than ac_acc. Breaks the inner loop. Default is 30.
#' @param ac_acc Numeric. Accuracy of the stopac break criterion in relation to the response. Default is 1/10000 of the function value at initial variables combination.
#' @export
#'
set_hyperparameters <- function(
    vf = NULL, rf = 1, dyn_rf = TRUE, t0 = 1e3, nlimit = 100,
    r = 0.6, k = 1, t_min = 0.1, maxgood = 100, stopac = 30, ac_acc = 1e-4
) {
  if (!is.null(vf) && !is.function(vf)) cli::cli_abort("vf must be a function or NULL")
  check_arg_len(1, rf, dyn_rf, t0, nlimit, r, k, t_min, maxgood, stopac, ac_acc)
  out <- named_list(vf, rf, dyn_rf, t0, nlimit, r, k, t_min, maxgood, stopac, ac_acc)
  return(out)
}

#' Boundary Constrainst for Simulated Annealing Optimization
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

#' Set Starting Values for Loss Function Optimization via Simulated Annealing
#' @description
#' Construct a named list of starting values for the optimization algorithm.
#' @param phi_n numeric scalar, specifying the starting value for the dispersion parameter of the negative binomial distribution
#' @param mu_n numeric scalar, specifying the starting value for the mean of the negative binomial distribution
#' @param mu_d numeric scalar, specifying the starting value for the mean of the normal distribution
#' @param sigma2_d numeric scalar, specifying the starting value for the standard deviation of the normal distribution
#' @param delta_hat numeric scalar, specifying the starting value for the expected effect size for sample size planning
#' @param w_pbs numeric scalar (between 0 and 1), specifying the starting value for the publiation bias weights
#' @return A named list of starting values for the optimization algorithm.
#' @export
set_start <- function(phi_n = NULL, mu_n = NULL, mu_d = NULL, sigma2_d = NULL, delta_hat = NULL, w_pbs = .5) {
  check_arg_len(0:1, phi_n, mu_n, mu_d, sigma2_d, delta_hat, w_pbs)
  named_list(phi_n, mu_n, mu_d, sigma2_d, delta_hat, w_pbs)
}

#' @title Estimate and Correct Publication Bias in Meta-Analysis under Consideration of Sample Size Planning
#' @param emp_data A dataframe or matrix with two columns: effect size (d) and sample size (n)
#' @param simetabias_control Control parameters \code{\link{simetabias_control}} for details.
#' @export
simetabias <- function(emp_data, simetabias_control = simetabias_control()) {
  # extract control parameters
  k_sim <- simetabias_control$k_sim
  bw <- simetabias_control$bw
  n_grid <- simetabias_control$n_grid
  alpha <- simetabias_control$alpha
  beta <- simetabias_control$beta
  slope_ssp <- simetabias_control$slope_ssp
  only_pb <- simetabias_control$only_pb
  trace <- simetabias_control$trace
  pr <- simetabias_control$pr
  # lists
  bounds <- simetabias_control$bounds
  start <- simetabias_control$start
  hyperparameters <- simetabias_control$hyperparameters
  only_pbs <- simetabias_control$only_pbs

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

  # which args not provided
  not_provided_start <- names(which(sapply(start, is.null)))
  if (any(not_provided_start %in% c("phi_n", "mu_n"))) {
    start_nb <- estimate_nb(emp_data[, "n"])
    if (any(not_provided_start == "phi_n")) start$phi_n <- start_nb["phi_n"]
    if (any(not_provided_start == "mu_n")) start$mu_n <- start_nb["mu_n"]
  }

  if (any(not_provided_start %in% c("mu_d", "sigma2_d"))) {
    start_norm <- estimate_norm(emp_data[, "d"])
    if (any(not_provided_start == "mu_d"))  start$mu_d <- start_norm["mu_d"]
    if (any(not_provided_start == "sigma2_d")) start$sigma2_d <- start_norm["sigma2_d"]
  }

  if(any(not_provided_start %in% "w_pbs")) {
    start$w_pbs <- 0.5
  }

  if(!only_pbs) {
    if(any(not_provided_start %in% "delta_hat")) {
      m <- mean(emp_data[, "d"])
      if(m == 0) m <- 0.1
      start$delta_hat <- m
    }
  }

  start_vec <- unlist(start)

  # run optimization
  hyperparameters$ac_acc <- hyperparameters$ac_acc * loss_function(start_vec)
  opt <- optimization::optim_sa(
    fun = loss_function, start = start_vec,
    lower = bounds$lower, upper = bounds$upper,
    maximization = FALSE,
    trace = trace,
    control = hyperparameters
  )
  return(opt)
}
