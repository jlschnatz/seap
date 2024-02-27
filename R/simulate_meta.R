#' Simulate Meta-Analytical Data with Selection Bias (Publication Bias and Sample Size Planning)
#' @description
#' This function simulates primary study data for a meta-analysis, incorporating selection biases such as publication bias
#' and bias due to sample size planning The simulation is based on a population of potential primary studies with assumed
#' marginal distributions of effect sizes (normal distribution) and sample sizes (negative-binomial distribution).
#'
#' The selection bias is introduced through two mechanisms:
#' - Publication bias in introduced via a selection mechanism by which studies with p-values above a specified alpha-level threshold are selected with probability `w_pbs`.
#' - Sample size planning bias is introduced by a selection mechanism by which the marginal distribution of sample sizes is weighted with a sigmoidal function that depends on the expected true effect size `delta_hat` and the type II error rate `beta`.
#'
#' @param k_sim number of simulated studies (approximate due to selection bias)
#' @param phi_n size parameter of negative binomial distribution
#' @param mu_n mean parameter of negative binomial sample size distribution
#' @param mu_d mean parameter of normal effect size distribution
#' @param sigma2_d variance parameter of normal effect size distribution
#' @param delta_hat expected true effect size
#' @param w_pbs weight for probability of selection given p(select | p > alpha)
#' @param slope_ssp discrimination parameter for sigmoid function
#' @param beta type II error rate
#' @param alpha type I error rate
#' @param only_pbs logical, if TRUE, only publication bias is applied
#' @return matrix of simulated data with columns named n and d and row size depends on selection mechanism
#' @examples
#' # Not run
#' x <- simulate_meta(
#'   k_sim = 100, phi_n = 5, mu_n = 100, mu_d = 0, sigma2_d = 1,
#'   delta_hat = 0.2, w_pbs = 0.5, slope_ssp = 4,
#'   beta = 0.2, alpha = 0.05
#'   )
#' head(x)
#' plot(x)
#' @export
#'
simulate_meta <- function(
    k_sim, phi_n, mu_n, mu_d, sigma2_d, delta_hat,
    w_pbs, slope_ssp = 4, beta = 0.2, alpha = 0.05, only_pbs = FALSE
) {
  if (any(c(phi_n, mu_n, sigma2_d) <= 0)) cli::cli_abort("phi_n, mu_n and sigma2_d must be positive")
  if (delta_hat == 0) cli::cli_abort("{.arg delta_hat} must be non-zero")
  if (beta < 0 | beta > 1) cli::cli_abort("{.arg beta} must be between 0 and 1")
  if (alpha < 0 | alpha > 1) cli::cli_abort("{.arg alpha} must be between 0 and 1")
  # first round
  first <- simulate_meta_(k_sim, phi_n, mu_n, mu_d, sigma2_d, delta_hat, w_pbs, slope_ssp, beta, alpha, only_pbs)
  # get proportion of survivors
  prop_surv <- sum(first$selected) / k_sim
  k_sim_adj <- ceiling(k_sim / prop_surv)
  if (is.infinite(k_sim_adj)) {
    matdist <- cbind(first$n, first$d)
    colnames(matdist) <- c("n", "d")
    return(matdist[first$selected, ])
  } else {
    # second round
    second <- simulate_meta_(k_sim_adj, phi_n, mu_n, mu_d, sigma2_d, delta_hat, w_pbs, slope_ssp, beta, alpha, only_pbs)
    # apply bias
    matdist <- cbind(second$n, second$d)
    colnames(matdist) <- c("n", "d")
    return(matdist[second$selected, ])
  }
}

simulate_meta_ <- function(k_sim, phi_n, mu_n, mu_d, sigma2_d, delta_hat,
                           w_pbs, slope_ssp = 4, beta = 0.2, alpha = 0.05, only_pbs = FALSE) {
  n <- stats::rnbinom(k_sim, size = phi_n, mu = mu_n)
  n[n == 0] <- 1
  gamma <- adjust_var(n, sigma2_d)
  d <- stats::rnorm(n = k_sim, mean = mu_d, sd = sqrt(gamma))
  p <- p_from_nd(n, d)
  bool_pbs <- pbs_select(p, w_pbs, alpha)
  if (only_pbs) {
    return(named_list(n, d, selected = bool_pbs))
  } else {
    bool_ssa <- ssa_select(n, delta_hat, slope_ssp, beta)
    return(named_list(n, d, selected = bool_pbs & bool_ssa))
  }
}
