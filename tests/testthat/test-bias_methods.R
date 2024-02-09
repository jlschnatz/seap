# Publication bias
test_that("Absent Publication Bias works", {
  expect_equal(pbs_select(p = 0.1, w_pbs = 1, alpha = 0.05), TRUE)
})

test_that("Complete Publication Bias works", {
  expect_equal(pbs_select(p = 0.1, w_pbs = 0, alpha = 0.05), FALSE)
})

test_that("Probabilistic Publication Bias works", {
  expect_equal(
    object = mean(pbs_select(p = rep(0.06, 1e6), w_pbs = 0.5, alpha = 0.05)),
    expected = 0.5,
    tolerance = 0.01
  )
})

# Sample size adjustment bias

test_that("Error if expected delta is 0", {
  expect_error(
    object = mean(ssa_select(n = rep(20, 1e6), d_delta_hat = 0, beta = 0.2)),
  )
})



