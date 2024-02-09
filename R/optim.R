# optimize_sa <- function(
#     start, compare_control = compare_control(), control_sa = list(),
#     maximization = FALSE, trace = FALSE, lower, upper
# ) {
#
#   fun <- function(params) compare_dists(params, control = compare_control)
#   opt <- optimization::optim_sa(
#     fun = fun, start = start, control = control_sa, maximization = maximization,
#     trace = trace, lower = lower, upper = upper
#   )
#   return(opt)
# }
#
#
#
# emp_dist <- simulate_meta(
#   k = 1000, k_size = 10, k_mu = 200, delta = 0.3, sigma2 = 0.2, d_delta_hat = 0.5, w_pbs = 1
# )
#
#
#
#
#
# t1 <- Sys.time()
# optimize_sa(
#   start = c(5, 150, 0.2, 0.3, 0.2, 0.5),
#   compare_control = compare_control(emp_dist = emp_dist, k = 1e4, n_rep = 10, parallel = TRUE, workers = 4),
#   control_sa = list(t0 = 100,
#                     nlimit = 100,
#                     t_min = 0.1,
#                     dyn_rf = FALSE,
#                     rf = 1,
#                     r = 0.8),
#   trace = TRUE,
#   lower = c(1, 20, -4, 0.01, 0.05, 0),
#   upper = c(100, 500, 4, 1, 5, 1)
# )
# t2 <- Sys.time()
#
#
# difference <- t2 - t1
#
# .Last.value -> test
# plot(test)
