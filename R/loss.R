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
  # Case 1: Q -> 0, P -> 0, P/Q -> NaN
  id_PQ0 <- which(P == 0 & Q == 0)
  P[id_PQ0] <- .Machine$double.xmin
  Q[id_PQ0] <- .Machine$double.xmin
  quot <- P / Q
  # Case 2: P -> 0, log(0) -> -Inf
  quot[quot == 0] <- .Machine$double.xmin
  # Case 3: Q -> 0, P/Q -> Inf
  quot[is.infinite(quot)] <- .Machine$double.xmax
  return(sum(P * log(quot)))
}


