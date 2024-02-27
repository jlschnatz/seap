#' Discrete Weighting Function for Publication Bias
#' @param p Numeric vector of p-values
#' @param w_pbs Numeric scalar. Weight for probability of selection given P(p | p > alpha)
#' @param alpha Numeric scalar. Type I error rate
#' @param detailed Logical scalar. Return detailed output (default: FALSE)
#' @returns If detailed = FALSE, returns a vector of weights. If detailed = TRUE, returns a list with p, w, w_pbs, and alpha.
#'
discrete_weight <- function(p, w_pbs, alpha, detailed = FALSE) {
  w <- ifelse(p < alpha, 1, w_pbs)
  if (detailed) {
    return(named_list(p, w, w_pbs, alpha))
  } else {
    return(w)
  }
}

#' Sigmoid Weighting Function for Sample Size Adjustment
#' @param n vector of sample sizes
#' @param delta_hat expected effect size
#' @param slope_ssp discrimination parameter (defaults to 5)
#' @param alpha type I error rate (default = 0.05)
#' @param beta type II error rate (default = 0.2)
#' @param detailed return detailed output (default: FALSE)
#' @returns If detailed = FALSE, returns a vector of weights. If detailed = TRUE, returns a list with n, w_ssa, n_req, delta_hat, alpha, and beta.
#'
sigmoid_weight <- function(n, delta_hat, slope_ssp = 5, alpha = 0.05, beta = 0.2, detailed = FALSE) {
  n_req <- 2 * pwr::pwr.t.test(d = delta_hat, power = 1 - beta, sig.level = alpha)$n
  w <- 1 / (1 + exp(-slope_ssp / n_req * (n - n_req)))
  if (detailed) {
    return(list(n = n, w_ssa = w, n_req = n_req, delta_hat = delta_hat, alpha = alpha, beta = beta))
  } else {
    return(w)
  }
}

#' Selection Bias based on Publication Bias
#' @param p vector of p-values
#' @param w_pbs weights for probability of selection given p(select | p > alpha)
#' @param alpha type I error rate
#' @return boolean vector of selection
pbs_select <- function(p, w_pbs, alpha) {
  stopifnot(rlang::is_scalar_double(w_pbs), w_pbs <= 1 & w_pbs >= 0)
  stopifnot(rlang::is_scalar_double(alpha), alpha < 1 & alpha > 0)
  if (alpha > 0.05) cli::cli_alert_warning("Type I error rate is larger than 0.05. This is not recommended.")
  selection_prob <- discrete_weight(p, w_pbs, alpha, detailed = FALSE)
  p_unif <- stats::runif(length(p), 0, 1)
  select_after_psb <- selection_prob > p_unif
  select_after_psb[is.na(select_after_psb)] <- FALSE # p-values can be NA if df = 0 -> then no selection
  return(select_after_psb)
}

#' Selection Bias based on Sample Size Adjustment
#' @param n vector of sample sizes
#' @param delta_hat expected effect size
#' @param slope_ssp discrimination parameter
#' @param beta type II error rate
#' @return boolean vector of length n
ssa_select <- function(n, delta_hat, slope_ssp = 0.5, beta = 0.2) {
  stopifnot(is.numeric(n), all(n >= 0))
  stopifnot(rlang::is_scalar_double(delta_hat))
  stopifnot(rlang::is_scalar_double(beta), beta < 1 & beta > 0)
  if (beta > 0.2) cli::cli_alert_warning("Type II error rate is larger than 0.2. This is not recommended.")
  if (delta_hat == 0) cli::cli_abort("Expected effect size cannot be 0 for {.fn pwr::pwr.t.test}")
  selection_prob <- sigmoid_weight(n, delta_hat, slope_ssp = slope_ssp, beta = beta, detailed = FALSE)
  p_unif <- stats::runif(length(n), 0, 1)
  select_after_ssa <- selection_prob > p_unif
  return(select_after_ssa)
}
