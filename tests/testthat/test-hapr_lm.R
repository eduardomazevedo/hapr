test_that("hapr_lm estimates coefficients correctly", {
  set.seed(123)  # For reproducibility
  devtools::load_all()

  # Create fake data
  n <- 1e4

  var_v <- 1/3
  var_epsilon <- 1/3
  var_thetaw <- 1/3

  true_improvement_ratio <- 1 / (1 - var_epsilon)

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

  # Call the hapr_lm function
  fit <- hapr_lm(y, gc_normalized, w, true_improvement_ratio)

  # Extract beta coefficients
  beta <- fit$beta

  # Check expected values
  expect_true(abs(beta[['gc']] - 0.42) < 0.05, 
              sprintf("Beta[gc] = %.3f is not close to 0.42", beta[['gc']]))
  expect_true(abs(beta[['w1']] - 0.17) < 0.05,
              sprintf("Beta[w1] = %.3f is not close to 0.17", beta[['w1']]))
  
  # Check remaining coefficients are close to zero
  remaining_betas <- setdiff(names(beta), c("gc", "w1"))
  for (param in remaining_betas) {
    expect_true(abs(beta[[param]]) < 0.05, 
                sprintf("Beta[%s] = %.3f is not close to zero", param, beta[[param]]))
  }
})
