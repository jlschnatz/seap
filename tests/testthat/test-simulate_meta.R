test_that("simulate meta works", {
  expect_equal(
    object = nrow(simulate_meta(
      k = 10,
      k_size = 10,
      k_mu = 30,
      delta = 0,
      sigma2 = 0.001,
      d_delta_hat = 0.2,
      w_pbs = 0,
      discr = 0,
      beta = 0.2,
      alpha = 0.05
    )),
    expected = 0
  )
})

test_that("Throw error if expected effect size for power analysis is zero", {
  expect_error(
    object = simulate_meta(
      k = 100,
      k_size = 10,
      k_mu = 30,
      delta = 0.3,
      d_delta_hat = 0,
      sigma2 = 0.4,
      w_pbs = 0,
      discr = 5,
      beta = 0.2,
      alpha = 0.05
    )
  )
})






