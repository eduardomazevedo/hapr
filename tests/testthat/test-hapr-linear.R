library(testthat)
library(tidyverse)
library(hapr)

# Helper: Simulate mock dataset with configurable parameters
simulate_mock_dataset <- function(n,
                                  var_v = 1/3,
                                  var_epsilon = 1/3,
                                  var_thetaw = 1/3,
                                  beta_gf = 0.42,
                                  beta_w1 = 0.17) {
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
  
  # Generate continuous outcome
  y <- beta_gf * gf + rnorm(n) + beta_w1 * w$w1
  
  list(
    w = w,
    gc = gc_normalized,
    y = y,
    beta_w1 = beta_w1
  )
}

# ---- STRUCTURE AND TYPE CHECKS (use small n for speed) ----
test_that("hapr structure and types are correct (linear, fast)", {
  data <- simulate_mock_dataset(n = 100)
  
  # Fit the model
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "lm",
    improvement_ratio = 1 / (1 - 0.2)  # consistent with a testable var_epsilon
  )
  
  # Extract estimated beta coefficients
  beta_hat <- fit$coefficients$beta
  
  # Check that beta hat is a double and named
  expect_type(beta_hat, "double")
  expect_named(beta_hat)
  
  # Simulate new data from the model
  sim_data <- hapr_simulate(fit, w = data$w)
  expect_s3_class(sim_data, "data.frame")
  
  # Validate core columns exist and have expected types
  for (col in c("gf", "gc")) {
    expect_type(sim_data[[col]], "double")
    expect_false(all(is.na(sim_data[[col]])))
    expect_gt(sd(sim_data[[col]], na.rm = TRUE), 0)
  }
  
  # Check that covariates are preserved in the simulated data
  expect_true(all(c("w1", "w2") %in% names(sim_data)))
  
  # Predict from simulated data
  preds <- predict(fit, newdata = sim_data)
  
  # Ensure predictions are returned and valid
  expect_true(all(c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w") %in% names(preds)))
  for (col in c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w")) {
    expect_false(any(is.na(preds[[col]])))
    expect_type(preds[[col]], "double")
  }
})

# ---- NUMERICAL ACCURACY TESTS (use large n) ----
test_that("hapr estimates coefficients correctly across parameter grid (linear)", {
  param_grid <- expand.grid(
    n = c(10000, 10000),
    var_v = c(0.2, 0.3),
    var_epsilon = c(0.2, 0.3),
    var_thetaw = c(0.3, 0.4),
    beta_w1 = c(0.2, 0.17),
    beta_gf = c(0.5, 0.4),
    KEEP.OUT.ATTRS = FALSE
  )
  
  for (i in seq_len(nrow(param_grid))) {
    params <- param_grid[i, ]
    label <- paste("n =", params$n,
                   "var_v =", params$var_v,
                   "var_epsilon =", params$var_epsilon,
                   "var_thetaw =", params$var_thetaw,
                   "beta_w1 =", params$beta_w1)
    
    # Generate simulated dataset for given parameter settings
    data <- simulate_mock_dataset(
      n = params$n,
      var_v = params$var_v,
      var_epsilon = params$var_epsilon,
      var_thetaw = params$var_thetaw,
      beta_w1 = params$beta_w1,
      beta_gf = params$beta_gf
    )
    
    # Fit the model
    fit <- hapr(
      y = data$y,
      gc = data$gc,
      w = data$w,
      model_type = "lm",
      improvement_ratio = 1 / (1 - params$var_epsilon)
    )
    
    # Extract estimated beta coefficients
    beta_hat <- fit$coefficients$beta
    
    # Check coefficient accuracy
    err <- abs(beta_hat["w1"] - data$beta_w1)
    expect_lt(err, 0.07)
    
    # Check factor coefficients are close to 0
    expect_true("w2B" %in% names(beta_hat))
    expect_true("w2C" %in% names(beta_hat))
    expect_lt(abs(beta_hat["w2B"]), 0.05)
    expect_lt(abs(beta_hat["w2C"]), 0.05)
  }
})

# ---- SNAPSHOT TEST (use small n) ----
test_that("hapr print output is stable (linear)", {
  data <- simulate_mock_dataset(n = 100)
  
  # Fit the model
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "lm",
    improvement_ratio = 1.5
  )
  
  # Capture stable printed output
  expect_snapshot(print(fit))
})