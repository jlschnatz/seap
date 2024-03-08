
if (interactive()) {
  emp_data <- simulate_meta(
    k_sim = 1e3, phi_n = 5, mu_n = 150, mu_d = 0, sigma2_d = 0.25,
    delta_hat = 5, w_pbs = .1,
    slope_ssp  = 0, only_pbs = TRUE
  )

  control <- speec_control(
    bw = "sheather-jones",
    n_grid = c(2^7+1, 2^7+1),
    pr = c(0.005, 0.995),
    k_sim = 1e4,
    bounds = set_boundaries(
      phi_n = c(1, 20),
      mu_n = c(50, 400),
      mu_d = c(-2, 2),
      sigma2_d = c(0.01, 1),
      delta_hat = c(0.05, 3),
      w_pbs = c(0, 1)
    ),
    start = set_start(
      phi_n = NULL,
      mu_n = NULL,
      mu_d = NULL,
      sigma2_d = NULL,
      delta_hat = NULL,
      w_pbs = 0.5
    ),
    alpha = 0.05,
    beta = 0.2,
    slope_ssp = 4,
    only_pbs = TRUE,
    trace = TRUE,
    hyperparameters = set_hyperparameters(
      ac_acc = 1e-4,
      nlimit = 10,
      r = .9,
      maxgood = 100,
      t0 = 1e4,
      dyn_rf = TRUE,
      vf = NULL,
      rf = 1,
      k = 1,
      t_min = 0.1,
      stopac = 30
    )
  )

  t1 <- Sys.time()
  result <- speec(
    emp_data = emp_data,
    speec_control = control
  )
  t2 <- Sys.time()
  t2 - t1
  print(result)
}
