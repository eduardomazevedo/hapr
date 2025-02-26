library(testthat)
set.seed(123)  # For reproducibility
devtools::load_all()

test_that("hapr_probit returns expected coefficients", {
  # Create fake data
  n <- 1e4
  
  var_v <- 1/3
  var_epsilon <- 1/3
  var_thetaw <- 1/3
  
  w <- data.frame(
    w1 = rnorm(n),
    w2 = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
  
  v <- rnorm(n) * sqrt(var_v)
  epsilon <- rnorm(n) * sqrt(var_epsilon)
  
  gf <- w$w1 * sqrt(var_thetaw) + v
  gc <- gf + epsilon
  gc_normalized <- scale(gc) |> as.numeric()
  
  y <- 0.42 * gf + rnorm(n) + 0.17 * w$w1
  
  y_binary <- as.numeric(y > mean(y)) |> as.factor()
  
  # Call the hapr_probit function
  full_fit <- hapr_probit(y_binary, gc_normalized, w, improvement_ratio = 1.5)
  
  # Test that coefficients are as expected
  expect_equal(full_fit$beta[['gf']], 0.42, tolerance = 0.1)
  expect_equal(full_fit$beta[['w1']], 0.17, tolerance = 0.1)
  expect_equal(full_fit$beta[['w2B']], 0, tolerance = 0.1)
})
