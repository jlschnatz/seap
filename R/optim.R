optimize_sa <- function(
    start, compare_control = compare_control(), control_sa = list(),
    maximization = FALSE, trace = FALSE, lower, upper
) {

  fun <- function(params) compare_dists(params, control = compare_control)
  opt <- optimization::optim_sa(
    fun = fun, start = start, control = control_sa, maximization = maximization,
    trace = trace, lower = lower, upper = upper
  )
  return(opt)
}


