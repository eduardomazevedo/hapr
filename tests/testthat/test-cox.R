library(testthat)
library(survival)

test_that("Cox model estimates coefficients correctly", {
  set.seed(123)  # For reproducibility
  
  # Number of observations
  n <- 1e4  
  var_v <- 1/3
  var_epsilon <- 2/3

  true_improvement_ratio <- 1 / (1 - var_epsilon)

  # Simulating gf ~ N(0, 1/3)
  gf <- rnorm(n, mean = 0, sd = sqrt(var_v))

  # Simulating w ~ N(0, 1)
  w <- rnorm(n, mean = 0, sd = 1)

  # Simulating epsilon ~ N(0, 2/3)
  epsilon <- rnorm(n, mean = 0, sd = sqrt(var_epsilon))

  # Computing gc = gf + epsilon
  gc <- gf + epsilon

  # Hazard rate: exp(0.42 * gf + 0.17 * w) / 10
  hazard_rate <- exp(0.42 * gf + 0.17 * w) / 10

  # Survival time t ~ Exp(hazard_rate)
  t <- rexp(n, rate = hazard_rate)

  # Storing data in a dataframe
  sim_data <- data.frame(gf, w, gc, hazard_rate, t)

  # Fit Cox model
  fit <- coxph(Surv(t) ~ gc + gf + w, data = sim_data)
  
  # Extract estimated coefficients
  beta <- coef(fit)

  # Test that estimates are close to the expected values
  expect_equal(beta[["gf"]], 0.42, tolerance = 0.05)
  expect_equal(beta[["w"]], 0.17, tolerance = 0.05)
})
