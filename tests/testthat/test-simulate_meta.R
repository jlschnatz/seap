test_that("simulate meta works", {
  expect_equal(
    object = nrow(simulate_meta(
      k_sim  = 10,
      phi_n = 10,
      mu_n = 30,
      mu_d = 0,
      sigma2_d = 0.001,
      delta_hat = 0.2,
      w_pbs = 0,
      slope_ssp = 0,
      beta = 0.2,
      alpha = 0.05,
      only_pbs = FALSE
    )),
    expected = 0
  )
})

test_that("Throw error if expected effect size for power analysis is zero", {
  expect_error(
    object = simulate_meta(
      k_sim  = 100,
      phi_n  = 10,
      mu_n  = 30,
      mu_d  = 0.3,
      sigma2_d = 0.4,
      delta_hat  = 0,
      w_pbs = 0,
      slope_ssp = 5,
      beta = 0.2,
      alpha = 0.05,
      only_pbs = FALSE
    )
  )
})
