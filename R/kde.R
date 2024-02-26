#' Estimate the parameters of a alternative parametrized negative binomial distribution
#' @param x a numeric vector
#' @param eps a numeric scalar that determines the tolerance for convergence
#' @return a named numeric vector with the estimated parameters (phi, mu)
#' @description
#' Uses the method of moments to estimate the parameters of a negative binomial distribution
#' @examples
#' # Not run
#' x <- rnbinom(100, size = 5, mu = 100)
#' estimate_nb(x, eps = 1e-5)
#'
estimate_nb <- function(x, eps = 1e-4) {
  n <- length(x)
  mu <- mean(x)
  dfr <- n - 1
  phi <- n/sum((x/mu - 1)^2)
  del <- 1
  while (abs(del) > eps) {
    phi <- abs(phi)
    del <- (sum(((x - mu)^2/(mu + mu^2/phi))) - dfr)/sum((x - mu)^2/(mu + phi)^2)
    phi <- phi - del
  }
  return(c("phi" = phi, "mu" = mu))
}

#' Estimate the parameters of a normal distribution
#' @param x a numeric vector
#' @return a named numeric vector with the estimated parameters (mu, sigma)
#' @description
#' Return the maximum likelihood estimates of the parameters of a normal distribution
#' @examples
#' # Not run
#' x <- rnorm(100, mean = 0, sd = 0.5)
#' estimate_norm(x)
estimate_norm <- function(x) c("mu" = mean(x), "sigma" = sd(x))

quantile_nb <- function(quantile = c(0.005, 0.995), phi, mu) {
  lb <- stats::qnbinom(quantile[1], size = phi, mu = mu, lower.tail = TRUE)
  ub <- stats::qnbinom(quantile[2], size = phi, mu = mu, lower.tail = TRUE)
  q <- c(lb, ub)
  names(q) <- c("lb", "ub")
  return(q)
}

quantile_norm <- function(quantile = c(0.005, 0.995), mu, sigma) {
  lb <- stats::qnorm(quantile[1], mean = mu, sd = sigma, lower.tail = TRUE)
  ub <- stats::qnorm(quantile[2], mean = mu, sd = sigma, lower.tail = TRUE)
  q <- c(lb, ub)
  names(q) <- c("lb", "ub")
  return(q)
}

find_kde_limits <- function(x, quantile = c(0.005, 0.995)) {
  if (!is.matrix(x)) cli::cli_abort("x must be a matrix")
  if (!all(c("d", "n") %in% colnames(x))) cli::cli_abort("x must have columns with names d and n")
  par_norm <- estimate_norm(x[, "d"])
  par_nb <- estimate_nb(x[, "n"])
  q_norm <- quantile_norm(quantile = quantile, mu = par_norm[1], sigma = par_norm[2])
  q_nb <- quantile_nb(quantile = quantile, phi = par_nb[1], mu = par_nb[2])
  lim_d <- c(min(c(x[, "d"], q_norm)), max(c(x[, "d"], q_norm)))
  names(lim_d) <- c("lb_d", "ub_d")
  lim_n <- c(min(c(x[, "n"], q_nb)), max(c(x[, "n"], q_nb)))
  names(lim_n) <- c("lb_n", "ub_n")
  return(list(n = lim_n, d = lim_d))
}


#' Estimate a bivariate kernel density estimate
#' @param x a numeric matrix with two columns (n, d)
#' @param bw a character string specifying the bandwidth selection method
#' @param n_grid a integer vector (2L) specifying the number of equally spaced gridpoints in each direction (should be a power of 2 + 1)
#' @return a matrix with the estimated densities
#' @examples
#' # Not run
#' k <- 100
#' X <- cbind(d = rnorm(k, mean = 0, sd = 0.5), n = rnbinom(k, size = 5, mu = 100))
#' bivariate_kde(X, bw = "silverman", n_grid = c(2^8 + 1, 2^8 + 1))
bivariate_kde <- function(x, bw = c("silverman", "scott", "sheather-jones", "ucv", "bcv"), n_grid, lims) {
  if(missing(bw)) cli::cli_abort("You must specify a bandwidth selection method h")
  if(missing(x)) cli::cli_abort("You must specify the data matrix X")
  if(dim(x)[2L] != 2L) cli::cli_abort("X must have two columns")
  if(missing(n_grid)) cli::cli_abort("You must specify the number equally spaced gridpoints in each direction")
  if(length(n_grid) != 2L) cli::cli_abort("n_grid must be a integer vector of length 2")
  if(!all(pow2_add1(n_grid))) {
    cli::cli_alert_warning("n_grid should be a power of 2 plus 1")
    n_grid[!pow2_add1(n_grid)] <- 2^ceiling(log2(n_grid[!pow2_add1(n_grid)])) + 1
    cli::cli_alert_info("Choosing the next larger power of 2 plus 1 as {n_grid}")
  }
  rlang::arg_match(bw)
  h <- switch(bw,
              "silverman" = apply(x, 2, bw.nrd0),
              "scott" = apply(x, 2, bw.nrd),
              "sheather-jones" = apply(x, 2, bw.SJ),
              "ucv" = apply(x, 2, bw.ucv),
              "bcv" = apply(x, 2, bw.bcv))
  est_ks <- KernSmooth::bkde2D(
    x = x[, names(lims)],
    bandwidth=h[names(lims)],
    gridsize = n_grid,
    range.x = lims
  )
  return(est_ks$fhat)
}

