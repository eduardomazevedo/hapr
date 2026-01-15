library(testthat)
library(tidyverse)
library(hapr)

# Helper: Simulate mock dataset with configurable parameters
simulate_mock_dataset <- function(n,
                                  var_v = 1/3,
                                  var_epsilon = 1/3,
                                  var_thetaw = 1/3,
                                  beta_gf = 0.42,
                                  beta_w1 = 0.17,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  w <- data.frame(
    w1 = rnorm(n),
    w2 = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )

  v <- rnorm(n) * sqrt(var_v)
  epsilon <- rnorm(n) * sqrt(var_epsilon)
  gf <- w$w1 * sqrt(var_thetaw) + v
  gc <- gf + epsilon
  gc_normalized <- scale(gc) |> as.numeric()
  y <- beta_gf * gf + rnorm(n) + beta_w1 * w$w1

  list(
    y = y,
    gc = gc_normalized,
    w = w,
    gf = gf,
    true_improvement_ratio = 1 / (1 - var_epsilon),
    beta_w1 = beta_w1,
    beta_gf = beta_gf
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
test_that("hapr estimates all coefficients correctly across variance scenarios (linear)", {
  # Test three variance scenarios: var_epsilon = 1/3, 3/4, 9/10
  # Other variances set so total = 1
  variance_scenarios <- list(
    list(var_epsilon = 1/3, var_v = 1/3, var_thetaw = 1/3),
    list(var_epsilon = 3/4, var_v = 1/8, var_thetaw = 1/8),
    list(var_epsilon = 9/10, var_v = 1/20, var_thetaw = 1/20)
  )
  
  n <- 1e4
  beta_gf_true <- 0.42
  beta_w1_true <- 0.17
  tolerance <- 0.01
  
  for (scenario in variance_scenarios) {
    var_epsilon <- scenario$var_epsilon
    var_v <- scenario$var_v
    var_thetaw <- scenario$var_thetaw
    
    label <- paste("var_epsilon =", var_epsilon, 
                   "var_v =", var_v, 
                   "var_thetaw =", var_thetaw)
    
    # Generate simulated dataset
    data <- simulate_mock_dataset(
      n = n,
      var_v = var_v,
      var_epsilon = var_epsilon,
      var_thetaw = var_thetaw,
      beta_w1 = beta_w1_true,
      beta_gf = beta_gf_true
    )

    # Fit the model
    fit <- hapr(
      y = data$y,
      gc = data$gc,
      w = data$w,
      model_type = "lm",
      improvement_ratio = 1 / (1 - var_epsilon)
    )

    # Extract estimated beta coefficients
    beta_hat <- fit$coefficients$beta

    # Check all coefficient accuracy
    expect_true("gf" %in% names(beta_hat), 
                info = paste("Missing 'gf' coefficient in", label))
    expect_true("w1" %in% names(beta_hat), 
                info = paste("Missing 'w1' coefficient in", label))
    
    # For beta_gf, compare to oracle regression with normalized Gf
    # (HAPR estimates beta for normalized Gf, not unnormalized)
    gf_normalized <- scale(data$gf) |> as.numeric()
    oracle_fit <- lm(data$y ~ gf_normalized + data$w$w1)
    beta_gf_oracle <- coef(oracle_fit)["gf_normalized"]
    
    # Build comparison table for error messages
    comparison_table <- tibble(
      Parameter = c("gf", "w1", "w2B", "w2C"),
      True = c(beta_gf_oracle, beta_w1_true, 0, 0),
      Estimated = c(beta_hat[["gf"]], beta_hat[["w1"]], 
                     beta_hat[["w2B"]], beta_hat[["w2C"]]),
      Error = c(abs(beta_hat[["gf"]] - beta_gf_oracle),
                abs(beta_hat[["w1"]] - beta_w1_true),
                abs(beta_hat[["w2B"]]),
                abs(beta_hat[["w2C"]]))
    )
    
    # Check each coefficient with detailed error message
    err_gf <- abs(beta_hat[["gf"]] - beta_gf_oracle)
    if (err_gf >= tolerance) {
      fail(paste("\n=== Test Case:", label, "===\n",
                 "gf coefficient failed tolerance check:\n",
                 "  Error =", err_gf, "(tolerance =", tolerance, ")\n",
                 "  HAPR estimate =", beta_hat[["gf"]], "\n",
                 "  Oracle (normalized) =", beta_gf_oracle, "\n\n",
                 "Full comparison table:\n",
                 capture.output(print(comparison_table)) |> paste(collapse = "\n")))
    }
    
    err_w1 <- abs(beta_hat[["w1"]] - beta_w1_true)
    if (err_w1 >= tolerance) {
      fail(paste("\n=== Test Case:", label, "===\n",
                 "w1 coefficient failed tolerance check:\n",
                 "  Error =", err_w1, "(tolerance =", tolerance, ")\n",
                 "  True =", beta_w1_true, "\n",
                 "  Estimated =", beta_hat[["w1"]], "\n\n",
                 "Full comparison table:\n",
                 capture.output(print(comparison_table)) |> paste(collapse = "\n")))
    }

    # Check factor coefficients are close to 0
    expect_true("w2B" %in% names(beta_hat), 
                info = paste("Missing 'w2B' coefficient in", label))
    expect_true("w2C" %in% names(beta_hat), 
                info = paste("Missing 'w2C' coefficient in", label))
    
    err_w2b <- abs(beta_hat[["w2B"]])
    if (err_w2b >= tolerance) {
      fail(paste("\n=== Test Case:", label, "===\n",
                 "w2B coefficient failed tolerance check:\n",
                 "  Error =", err_w2b, "(tolerance =", tolerance, ")\n",
                 "  Estimated =", beta_hat[["w2B"]], "\n\n",
                 "Full comparison table:\n",
                 capture.output(print(comparison_table)) |> paste(collapse = "\n")))
    }
    
    err_w2c <- abs(beta_hat[["w2C"]])
    if (err_w2c >= tolerance) {
      fail(paste("\n=== Test Case:", label, "===\n",
                 "w2C coefficient failed tolerance check:\n",
                 "  Error =", err_w2c, "(tolerance =", tolerance, ")\n",
                 "  Estimated =", beta_hat[["w2C"]], "\n\n",
                 "Full comparison table:\n",
                 capture.output(print(comparison_table)) |> paste(collapse = "\n")))
    }
  }
})

# ---- SNAPSHOT TEST (use small n) ----
test_that("hapr print output is stable (linear)", {
  set.seed(123)
  data <- simulate_mock_dataset(n = 1000)

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

# ---- CI COVERAGE TEST ----
test_that("hapr confidence intervals for beta achieve exact match coverage (linear)", {
  skip_if(Sys.getenv("HAPR_RUN_COVERAGE_TESTS") == "", 
          "Skipping coverage test. Set HAPR_RUN_COVERAGE_TESTS=1 to run.")
  
  # ---- Part 1: Configuration ----
  default_params <- list(
    n = 1e4,
    var_v = 1 / 3,
    var_epsilon = 1 / 3,
    var_thetaw = 1 / 3,
    beta_gf = 0.42,
    beta_w1 = 0.17,
    seed = 123
  )

  create_simulated_dataset <- function(params = list()) {
    p <- modifyList(default_params, params)
    set.seed(p$seed)

    n <- p$n
    var_v <- p$var_v
    var_epsilon <- p$var_epsilon
    var_thetaw <- p$var_thetaw
    beta_gf <- p$beta_gf
    beta_w1 <- p$beta_w1

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
    y <- beta_gf * gf + rnorm(n) + beta_w1 * w$w1

    list(
      y = y,
      gc = gc_normalized,
      w = w,
      true_improvement_ratio = true_improvement_ratio
    )
  }

  # ---- Part 2: True beta from simulation config ----
  sim_data <- create_simulated_dataset()
  fit <- hapr(sim_data$y, sim_data$gc, sim_data$w,
              model_type = "lm",
              improvement_ratio = sim_data$true_improvement_ratio)

  beta_names <- names(fit$coefficients$beta)

  # True parameter
  true_beta <- c(
    `(Intercept)` = 0,
    gf = default_params$beta_gf,
    w1 = default_params$beta_w1,
    w2B = 0,
    w2C = 0
  )

  n_sim <- 1000
  covered_matrix <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
  colnames(covered_matrix) <- beta_names

  # ---- Part 3: Run simulations ----
  for (i in 1:n_sim) {
    if (i %% 50 == 0) cat("Simulation", i, "\n")

    sim_data_i <- create_simulated_dataset(params = list(seed = i))
    fit_i <- hapr(
      y = sim_data_i$y,
      gc = sim_data_i$gc,
      w = sim_data_i$w,
      model_type = "lm",
      improvement_ratio = sim_data_i$true_improvement_ratio
    )

    ci_beta <- fit_i$ci_beta

    for (term in beta_names) {
      ci_lower <- ci_beta[term, "Lower"]
      ci_upper <- ci_beta[term, "Upper"]
      covered_matrix[i, term] <- (true_beta[term] >= ci_lower) & (true_beta[term] <= ci_upper)
    }
  }

  # ---- Part 4: Compute coverage ----
  coverage_df <- tibble(
    Term = beta_names,
    Coverage = round(colMeans(covered_matrix), 3)
  )

  print(coverage_df)

  expect_true(all(coverage_df$Coverage > 0.90),
              info = paste("Coverage below threshold:\n",
                           paste(coverage_df$Term, coverage_df$Coverage, collapse = "\n")))
})

