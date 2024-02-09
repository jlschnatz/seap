#' Discrete Weighting Function for Publication Bias
#' @param p vector of p-values
#' @param w_pbs weights for probability of selection given P(p | p > alpha)
#' @param alpha type I error rate
#' @param detailed return detailed output (default: FALSE)
#' @returns If detailed = FALSE, returns a vector of weights. If detailed = TRUE, returns a list with p, w, w_pbs, and alpha.
#' @export
discrete_weight <- function(p, w_pbs, alpha, detailed = FALSE) {
  w <- ifelse(p < alpha, 1, w_pbs)
  if (detailed) {
    return(list(p = p, w = w, w_pbs = w_pbs, alpha = alpha))
  } else {
    return(w)
  }
}

#' Sigmoid Weighting Function for Sample Size Adjustment
#' @param n vector of sample sizes
#' @param d_delta_hat expected effect size
#' @param discr discrimination parameter (defaults to 5)
#' @param alpha type I error rate (default = 0.05)
#' @param beta type II error rate (default = 0.2)
#' @param detailed return detailed output (default: FALSE)
#' @returns If detailed = FALSE, returns a vector of weights. If detailed = TRUE, returns a list with n, w_ssa, n_req, delta_hat, alpha, and beta.
#' @export
sigmoid_weight <- function(n, d_delta_hat, discr = 5, alpha = 0.05, beta = 0.2, detailed = FALSE) {
  n_req <- 2*pwr::pwr.t.test(d = d_delta_hat, power = 1 - beta, sig.level = alpha)$n
  w <- 1 / (1 + exp(-discr/n_req * (n - n_req)))
  if (detailed) {
    return(list(n = n, w_ssa = w, n_req = n_req, d_delta_hat = d_delta_hat, alpha = alpha, beta = beta))
  } else {
    return(w)
  }
}

#' Selection Bias based on Publication Bias
#' @param p vector of p-values
#' @param w_pbs weights for probability of selection given p(select | p > alpha)
#' @param alpha type I error rate
#' @return boolean vector of selection
#' @examples
#' p <- c(0.02, 0.06, 0.1)
#' pbs_select(p, 0.2, 0.05)
#' @export
#'
pbs_select <- function(p, w_pbs, alpha) {
  stopifnot(is_scalar(w_pbs), w_pbs <= 1 & w_pbs >= 0)
  stopifnot(is_scalar(alpha), alpha < 1 & alpha > 0)
  if (alpha > 0.05) cli::cli_alert_warning("Type I error rate is larger than 0.05. This is not recommended.")
  selection_prob <- discrete_weight(p, w_pbs, alpha, detailed = FALSE)
  p_unif <- stats::runif(length(p), 0, 1)
  select_after_psb <- selection_prob > p_unif
  select_after_psb[is.na(select_after_psb)] <- FALSE # p-values can be NA if df = 0 -> then no selection
  return(select_after_psb)
}

#' Selection Bias based on Sample Size Adjustment
#' @param n vector of sample sizes
#' @param d_delta_hat expected effect size
#' @param discr discrimination parameter
#' @param beta type II error rate
#' @return boolean vector of selection
#' @examples
#' n <- c(50, 20, 130)
#' ssa_select(n, 0.2, 0.2)
#' @export
#'
ssa_select <- function(n, d_delta_hat, discr = 0.5, beta = 0.2) {
  stopifnot(is.numeric(n), all(n >= 0))
  stopifnot(is_scalar(d_delta_hat))
  stopifnot(is_scalar(beta), beta < 1 & beta > 0)
  if (beta > 0.2) cli::cli_alert_warning("Type II error rate is larger than 0.2. This is not recommended.")
  if (d_delta_hat == 0) cli::cli_abort("Expected effect size cannot be 0 for `pwr::pwr.t.test()`")
  selection_prob <- sigmoid_weight(n, d_delta_hat, discr = discr, beta = beta, detailed = FALSE)
  p_unif <- stats::runif(length(n), 0, 1)
  select_after_ssa <- selection_prob > p_unif
  return(select_after_ssa)
}

