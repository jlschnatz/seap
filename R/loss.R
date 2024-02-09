#' Low-level KL divergence computation
#' @description
#' This function computes the KL-divergence for discrete distributions. It handles cases where due to numerical precision,
#' p -> 0 and log(0) -> -Inf and Q -> 0 and P/Q -> Inf. In these cases, the function replaces 0 with the smallest positive number
#' and Inf with the largest positive number.
#' @param P empirical distribution
#' @param Q theoretical distribution
#' @return numeric value of the KL-divergence
kl_div <- function(P, Q) {
  if (any(dim(P) != dim(Q))) cli::cli_abort("P and Q must have the same dimensions")
  quot <- P / Q
  # Case 1: P -> 0, log(0) -> -Inf
  quot[quot == 0] <- .Machine$double.xmin
  # Case 2: Q -> 0, P/Q -> Inf
  quot[is.infinite(quot)] <- .Machine$double.xmax
  return(sum(P * log(quot)))
}

#' Compute the Kullback-Leibler divergence between the empirical and theoretical effect size sample size joint distribution
#' @description
#' This function computes the Kullback-Leibler divergence between the empirical and theoretical effect size sample size joint distribution.
#' Empirical and theoretical distributions are estimated using 2-dimensional kernel density estimation with a grid size set by `n_grid`.
#' The lower and upper bounds for the grid are computed based on the absolute minimum and maximum of the data and the quantiles of the fitted distributions.
#' The starting values for the maximum likelihood estimation of the normal and negative binomial distributions can be set with `start_norm` and `start_nb`.
#' @param p data frame or matrix of the empirical distribution with colnames n and d
#' @param q data frame or matrix of the theoretical distribution with colnames n and d
#' @param n_grid integer with the number of grid points (k^2) for the 2-dimensional kernel density estimation process
#' @param bounds vector of length 2 with lower and upper quantiles (default: (0.005, 0.995))
#' @param start_nb vector of length 2 with starting values for ML-estimation for phi and mu of negative binomial distribution (default: (10, 100)
#' @param start_norm vector of length 2 with starting values for ML-estimation for mu and sigma of normal distribution (default: (0, 0.5)
#' @return numeric value of the Kullback-Leibler divergence
#' @export
compute_loss <- function(p, q, n_grid = 100, bounds = c(0.005, 0.995), start_nb = c(10, 100), start_norm = c(0, 0.5)) {
  if (!all(colnames(p) %in% c("n", "d") & colnames(p) == colnames(q))) cli::cli_abort("p and q must have the same column names (n, d)")
  lims <- find_kde_bounds(x = p, bounds = bounds, start_nb = start_nb, start_norm = start_norm)
  p_dens_mat <- MASS::kde2d(p[, "n"], p[, "d"], n = n_grid, lims = lims)$z
  w_dens_mat <- MASS::kde2d(p[, "n"], q[, "d"], n = n_grid, lims = lims)$z
  return(kl_div(p_dens_mat, w_dens_mat))
}

#' Grid Bounds for 2-dimensional KDE
#' @description
#' Computes the lower and upper bounds for the grid of the 2-dimensional KDE based on the absolute
#' minimum and maximum of the data and the quantiles of the fitted distributions.
#' @param x data frame or matrix with colnames d and n
#' @param bounds vector of length 2 with lower and upper quantiles
#' @param start_norm vector of length 2 with starting values for ML-estimation for mu and sigma of normal distribution
#' @param start_nb vector of length 2 with starting values for ML-estimation for phi and mu of negative binomial distribution
#' @return named vector of length 4 with lower and upper bounds for d and n
#' @export
#' @seealso [compute_loss()]
#' @examples
#' x <- cbind(d = rnorm(100, mean = 0.2, sd = 0.4), n = rnbinom(100, size = 10, mu = 100))
#' find_kde_bounds(x = x, bounds = c(0.005, 0.995), start_norm = c(0, 1), start_nb = c(10, 100))
find_kde_bounds <- function(x, bounds = c(0.005, 0.995), start_norm = c(0, 1), start_nb = c(10, 100)) {
  if (!is.matrix(x)) cli::cli_abort("x must be a matrix")
  if (!all(c("d", "n") %in% colnames(x))) cli::cli_abort("x must have columns with names d and n")
  opt_norm <- optim_norm(x = x[, "d"], mu = start_norm[1], sigma = start_norm[2])
  opt_nb <- optim_nb(x = x[, "n"], phi = start_nb[1], mu = start_nb[2])
  q_norm <- norm_quantile(bounds = bounds, mu = opt_norm[1], sigma = opt_norm[2])
  q_nb <- nb_quantile(bounds = bounds, phi = opt_nb[1], mu = opt_nb[2])
  lim_d <- c(min(c(x[, "d"], q_norm)), max(c(x[, "d"], q_norm)))
  names(lim_d) <- c("lb_d", "ub_d")
  lim_n <- c(min(c(x[, "n"], q_nb)), max(c(x[, "n"], q_nb)))
  names(lim_n) <- c("lb_n", "ub_n")
  return(c(lim_n, lim_d))
}

optim_norm <- function(x, mu, sigma) {
  nll <- function(pars, data) {
    mu <- pars[1]
    sigma <- pars[2]
    -sum(stats::dnorm(x = data, mean = mu, sd = sigma, log = TRUE))
  }
  opt <- suppressWarnings(stats::optim(par = c(mu, sigma), nll, data = x))
  names(opt$par) <- c("mu", "sigma")
  return(opt$par)
}

norm_quantile <- function(bounds = c(0.005, 0.995), mu, sigma) {
  lb <- stats::qnorm(bounds[1], mu, sigma, lower.tail = TRUE)
  ub <- stats::qnorm(bounds[2], mu, sigma, lower.tail = TRUE)
  q <- c(lb, ub)
  names(q) <- c("lb", "ub")
  return(q)
}

optim_nb <- function(x, phi, mu) {
  nll <- function(pars, data) {
    phi <- pars[1]
    mu <- pars[2]
    -sum(stats::dnbinom(x = data, mu = mu, size = phi, log = TRUE))
  }
  opt <- suppressWarnings(stats::optim(par = c(phi, mu), nll, data = x))
  names(opt$par) <- c("phi", "mu")
  return(opt$par)
}

nb_quantile <- function(bounds = c(0.005, 0.995), phi, mu) {
  lb <- stats::qnbinom(bounds[1], size = phi, mu = mu, lower.tail = TRUE)
  ub <- stats::qnbinom(bounds[2], size = phi, mu = mu, lower.tail = TRUE)
  q <- c(lb, ub)
  names(q) <- c("lb", "ub")
  return(q)
}
