#' Simulate Meta-Analytical Data with Selection Bias (Publication Bias and Sample Size Adjustment Bias)
#' @description
#' This function simulates primary study data for a meta-analysis, incorporating selection biases such as publication bias
#' and bias due to sample size adjustment. The simulation is based on a population of potential primary studies with assumed
#' marginal distributions of effect sizes (normal distribution) and sample sizes (negative-binomial distribution).
#'
#' The selection bias is introduced through two mechanisms:
#'
#' - Publication bias is introduced by a selection mechanism that chooses studies with p-values above a specified alpha-level
#'   with a probability denoted by `w_pbs`.
#'
#' - Sample size adjustment bias is introduced by a selection mechanism that chooses studies with sample sizes below a required
#'   sample size, given estimates from power-analysis with a specified power-level and expected effect size. The probability
#'   of selection is denoted by `w_ssa`.
#'
#' @param k number of primary studies in the 'population'
#' @param n_size size parameter of negative binomial distribution
#' @param n_mu mean parameter of negative binomial sample size distribution
#' @param d_mu mean parameter of normal effect size distribution
#' @param d_sigma variance parameter of normal effect size distribution
#' @param d_delta_hat expected true effect size
#' @param w_pbs weight for probability of selection given p(select | p > alpha)
#' @param discr discrimination parameter for sigmoid function
#' @param beta type II error rate
#' @param alpha type I error rate
#' @return matrix of simulated data with columns named n and d and row size depends on selection mechanism
#' @export
#'
simulate_meta <- function(k, n_size, n_mu, d_mu, d_sigma, d_delta_hat, w_pbs, discr = 5, beta = 0.2, alpha = 0.05) {
  requireNamespace("cli", quietly = TRUE)
  if (any(c(n_size, n_mu, d_sigma) <= 0)) cli::cli_abort("n_size, n_mu and d_sigma must be positive")
  if (d_delta_hat == 0)  cli::cli_abort("d_delta_hat must be non-zero")
  if (beta < 0 | beta > 1) cli::cli_abort("beta must be between 0 and 1")
  if (alpha < 0 | alpha > 1) cli::cli_abort("alpha must be between 0 and 1")
  # first round
  n1 <- stats::rnbinom(k, size = n_size, mu = n_mu)
  n1[n1 == 0] <- 1
  se1 <- rescale_se_d(n1, d_sigma)
  d1 <- stats::rnorm(n = k, mean = d_mu, sd = se1)
  p1 <- p_from_nd(n1, d1)
  bool_pbs1 <- pbs_select(p1, w_pbs, alpha)
  bool_ssa1 <- ssa_select(n1, d_delta_hat, discr, beta)
  # get proportion of survivors
  prop_surv <- sum(bool_pbs1 & bool_ssa1)/k
  k_adj <- ceiling(k/prop_surv)
  if (is.infinite(k_adj)) {
    matdist <- cbind(n1, d1)
    colnames(matdist) <- c("n", "d")
    return(matdist[bool_pbs1 & bool_ssa1, ])
  } else {
  # second round
  n2 <- stats::rnbinom(k_adj, size = n_size, mu = n_mu)
  n2[n2 == 0] <- 1
  se2 <- rescale_se_d(n2, d_sigma)
  d2 <- stats::rnorm(n = k_adj, mean = d_mu, sd = se2)
  p2 <- p_from_nd(n2, d2)
  bool_pbs2 <- pbs_select(p2, w_pbs, alpha)
  bool_ssa2 <- ssa_select(n2, d_delta_hat, discr, beta)
  # apply bias
  matdist <- cbind(n2, d2)
  colnames(matdist) <- c("n", "d")
  return(matdist[bool_pbs2 & bool_ssa2, ])
  }
}

