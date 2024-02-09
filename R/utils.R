#' Calculate p-value from sample size and effect size
#' @details
#' This calculation assumes a two-sided two-sample t-test with equal variance
#' @param n vector of sample sizes
#' @param d vector of effect sizes
#' @return vector of p-values
#'
p_from_nd <- function(n, d) {
  n_vec <- get_2n(n)
  t <- d / sqrt(1/n_vec[, 1] + 1/n_vec[, 2])
  df <- sum(n_vec) - 2
  p <- suppressWarnings(2 * stats::pt(abs(t), df, lower.tail = FALSE))
  return(p)
}

#' Calculate p-value from sample size and effect size
#' @details
#' This calculation assumes a two-sided two-sample t-test with equal variance
#' @param n vector of sample sizes
#' @param lower.tail logical vector of 1L, FALSE indicates positive effect sizes and TRUE indicates negative effect sizes
#' @param alpha type I error rate
#' @return vector of d-values
#'
d_from_np <- function(n, lower.tail = FALSE, alpha = 0.05) {
  n_vec <- get_2n(n)
  df <- sum(n_vec) - 2
  t_value <- stats::qt(alpha/2, df, lower.tail = lower.tail)
  d <- t_value * sqrt(1/n_vec[1] + 1/n_vec[2])
  return(d)
}

is_scalar <- function(x) is.numeric(x) && length(x) == 1

get_2n <- function(n) {
  stopifnot(is.numeric(n))
  cbind(ifelse(n %% 2 == 0, n/2, ceiling(n/2)), ifelse(n %% 2 == 0, n/2, floor(n/2)))
}

#' Selfnaming list
#' Taken from lme4::namedList()
#' @param ... comma-separated arguments
#' @usage named_list(...)
named_list <- function(...) {
  l <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1]
  nm <- names(l)
  if (is.null(nm)) nm <- snm
  nonames <- nm == ""
  if (any(nonames)) nm[nonames] <- snm[nonames]
  names(l) <- nm
  return(l)
}

