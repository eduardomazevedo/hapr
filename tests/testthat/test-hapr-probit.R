library(testthat)
library(tidyverse)
library(hapr)

# Helper: Simulate mock dataset
simulate_mock_dataset <- function(n = 1e5, var_v = 1/3, var_epsilon = 1/3, var_thetaw = 1/3) {
  set.seed(123)
  
  # Generate covariates
  w <- data.frame(
    w1 = rnorm(n),
    w2 = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )
  
  # Generate latent and observed variables
  v <- rnorm(n) * sqrt(var_v)
  epsilon <- rnorm(n) * sqrt(var_epsilon)
  gf <- w$w1 * sqrt(var_thetaw) + v
  gc <- gf + epsilon
  gc_normalized <- scale(gc) |> as.numeric()
  
  # Generate binary outcome
  y <- 0.42 * gf + rnorm(n) + 0.17 * w$w1
  y_binary <- as.numeric(y > mean(y)) |> as.factor()
  
  list(w = w, gc = gc_normalized, y_binary = y_binary)
}

# Main test
test_that("hapr estimates coefficients close to true values", {
  data <- simulate_mock_dataset()
  
  # Fit the model
  fit <- hapr(
    y = data$y_binary,
    gc = data$gc,
    w = data$w,
    model_type = "probit",
    improvement_ratio = 1.5
  )
  
  # Extract estimated beta coefficients
  beta_hat <- fit$coefficients$beta
  
  # Check that beta hat is a double
  expect_type(beta_hat, "double")
  expect_named(beta_hat)
  
  # We want to check that beta_hat["w1"] is very close to the real value 0.17.
  expect_true(abs(beta_hat["w1"] - 0.17) < 0.05)
  
  
  expect_true("w2B" %in% names(beta_hat))
  expect_true("w2C" %in% names(beta_hat))
  
  # Assert factor coefficients are near zero (not predictive in data)
  expect_true(abs(beta_hat["w2B"]) < 0.05)
  expect_true(abs(beta_hat["w2C"]) < 0.05)
  
  # Check that coefficients are not extreme
  expect_true(all(abs(beta_hat) < 10))
})


test_that("predictions from simulated data are internally consistent", {
  data <- simulate_mock_dataset()
  
  # Fit model
  fit <- hapr(
    y = data$y_binary,
    gc = data$gc,
    w = data$w,
    model_type = "probit",
    improvement_ratio = 1.5
  )
  
  # Simulate new data from the model
  sim_data <- hapr_simulate(fit, w = data$w)
  expect_s3_class(sim_data, "data.frame")
  required_cols <- c("gf", "gc")
  # Validating simulation columns
  for (col in required_cols) {
    expect_type(sim_data[[col]], "double")
    expect_false(all(is.na(sim_data[[col]])))
    expect_gt(sd(sim_data[[col]], na.rm = TRUE), 0)
  }
  
  # Check for covariates
  expect_true("w1" %in% names(sim_data))
  expect_true("w2" %in% names(sim_data))
  # Predict from simulated data
  preds <- predict(fit, newdata = sim_data)
  
  # Ensure predictions are valid probabilities
  for (col in c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w")) {
    expect_true(all(preds[[col]] >= 0 & preds[[col]] <= 1, na.rm = TRUE))
    expect_false(any(is.na(preds[[col]])))
  }
  # Ensure predictions exist
  expect_true(all(c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w") %in% names(preds)))
  
  # Internal consistency: predictions from w and gf should be correlated
  expect_gt(cor(preds$y_hat_w, preds$y_hat_gf_w, use = "complete.obs"), 0.8)
  
  # w and gc-based predictions should be correlated too
  expect_gt(cor(preds$y_hat_w, preds$y_hat_gc_w, use = "complete.obs"), 0.8)
  
  # gc and gf based predictions should be correlated too
  expect_gt(cor(preds$y_hat_gc_w, preds$y_hat_gf_w, use = "complete.obs"), 0.7)
  
  # y_hat_gf_w increases with gf
  cor_gf <- cor(preds$y_hat_gf_w, sim_data$gf, use = "complete.obs")
  expect_gt(cor_gf, 0.8)
})

test_that("hapr print output is stable", {
  data <- simulate_mock_dataset()

  fit <- hapr(
    y = data$y_binary,
    gc = data$gc,
    w = data$w,
    model_type = "probit",
    improvement_ratio = 1.5
  )

  expect_snapshot(print(fit))
})
