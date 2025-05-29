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
  
  # Continuous outcome for linear model
  y <- 0.42 * gf + rnorm(n) + 0.17 * w$w1
  
  list(w = w, gc = gc_normalized, y = y)
}

# Main test
test_that("hapr estimates coefficients close to true values (linear)", {
  data <- simulate_mock_dataset()
  
  # Fit the model
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "lm",
    improvement_ratio = 1.5
  )
  
  beta_hat <- fit$coefficients$beta
  
  expect_type(beta_hat, "double")
  expect_named(beta_hat)
  
  expect_true(abs(beta_hat["w1"] - 0.17) < 0.05)
  
  expect_true("w2B" %in% names(beta_hat))
  expect_true("w2C" %in% names(beta_hat))
  
  expect_true(abs(beta_hat["w2B"]) < 0.05)
  expect_true(abs(beta_hat["w2C"]) < 0.05)
  
  expect_true(all(abs(beta_hat) < 10))
})

test_that("predictions from simulated data are internally consistent (linear)", {
  data <- simulate_mock_dataset()
  
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "lm",
    improvement_ratio = 1.5
  )
  
  sim_data <- hapr_simulate(fit, w = data$w)
  expect_s3_class(sim_data, "data.frame")
  required_cols <- c("gf", "gc")
  
  for (col in required_cols) {
    expect_type(sim_data[[col]], "double")
    expect_false(all(is.na(sim_data[[col]])))
    expect_gt(sd(sim_data[[col]], na.rm = TRUE), 0)
  }
  
  expect_true("w1" %in% names(sim_data))
  expect_true("w2" %in% names(sim_data))
  
  preds <- predict(fit, newdata = sim_data)
  
  expect_true(all(c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w") %in% names(preds)))
  
  # Skip [0, 1] checks since linear predictions are not probabilities
  for (col in c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w")) {
    expect_false(any(is.na(preds[[col]])))
    expect_type(preds[[col]], "double")
  }
  
  expect_gt(cor(preds$y_hat_w, preds$y_hat_gf_w, use = "complete.obs"), 0.8)
  expect_gt(cor(preds$y_hat_w, preds$y_hat_gc_w, use = "complete.obs"), 0.8)
  expect_gt(cor(preds$y_hat_gc_w, preds$y_hat_gf_w, use = "complete.obs"), 0.7)
  
  cor_gf <- cor(preds$y_hat_gf_w, sim_data$gf, use = "complete.obs")
  expect_gt(cor_gf, 0.8)
})

test_that("hapr print output is stable (linear)", {
  data <- simulate_mock_dataset()
  
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "lm",
    improvement_ratio = 1.5
  )
  
  expect_snapshot(print(fit))
})