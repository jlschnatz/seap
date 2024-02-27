#' Calculate p-value from sample size and effect size
#' @details
#' This calculation assumes a two-sided two-sample t-test with equal variance
#' @param n vector of sample sizes
#' @param d vector of effect sizes
#' @return vector of p-values
#'
p_from_nd <- function(n, d) {
  n_vec <- get_2n(n)
  t <- d / sqrt(1 / n_vec[, 1] + 1 / n_vec[, 2])
  df <- rowSums(n_vec) - 2
  p <- suppressWarnings(2 * stats::pt(abs(t), df, lower.tail = FALSE))
  return(p)
}

#' Split vector of total sample sizes into a matrix samples sizes by group
#' @param n vector of total sample sizes
#' @return matrix with two columns and len n rows of sample sizes
get_2n <- function(n) {
  stopifnot(is.numeric(n))
  cbind(ifelse(n %% 2 == 0, n / 2, ceiling(n / 2)), ifelse(n %% 2 == 0, n / 2, floor(n / 2)))
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

#' Adjust variance for sample size
#' @param n vector of sample sizes
#' @param sigma2_d variance of the effect size
#' @return vector of adjusted variances
adjust_var <- function(n, sigma2_d) {
  var_mean <- sigma2_d / n
  return(var_mean / mean(var_mean) * sigma2_d)
}

#' Check if a number is a power of 2 plus 1
#' @return a logical vector
#' @param x a integer vector
pow2_add1 <- function(x) {
  xint <- vctrs::vec_cast(x, to = integer())
  bitwAnd(xint - 1, xint - 2) == 0
}

#' Check arguments lengths
#' @param len integer vector of expected lengths
#' @param ... list of arguments
check_arg_len <- function(len, ...) {
  arg <- named_list(...)
  lens <- lengths(arg)
  not_match <- which(!lens %in% len)
  if(!rlang::is_empty(not_match)) {
    text <- cli::cli_vec(paste0(len, "L"), style = list("vec-last" = " or "))
    cli::cli_abort("Argument {names(arg)[not_match]} is not length {text}")
  }
}
