library(testthat)
library(tidyverse)
library(hapr)

# Helper: Simulate mock dataset
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

  # Binary outcome for probit model
  latent <- beta_gf * gf + rnorm(n) + beta_w1 * w$w1
  y <- as.factor(as.numeric(latent > 0))

  list(
    w = w,
    gc = gc_normalized,
    y = y,
    beta_w1 = beta_w1
  )
}

# ---- STRUCTURE AND TYPE CHECKS (use small n for speed) ----
test_that("hapr structure and types are correct (probit, fast)", {
  var_epsilon <- 0.2
  data <- simulate_mock_dataset(n = 100, var_epsilon = var_epsilon)

  # Fit the model
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "probit",
    improvement_ratio = 1 / (1 - var_epsilon)
  )

  # Extract estimated beta coefficients
  beta_hat <- fit$coefficients$beta

  # Check that beta hat is a double
  expect_type(beta_hat, "double")
  expect_named(beta_hat)

  # Simulate new data from the model
  sim_data <- hapr_simulate(fit, w = data$w)
  expect_s3_class(sim_data, "data.frame")

  # Validating simulation columns
  for (col in c("gf", "gc")) {
    expect_type(sim_data[[col]], "double")
    expect_false(all(is.na(sim_data[[col]])))
    expect_gt(sd(sim_data[[col]], na.rm = TRUE), 0)
  }

  # Check for covariates
  expect_true(all(c("w1", "w2") %in% names(sim_data)))

  # Predict from simulated data
  preds <- predict(fit, newdata = sim_data)

  # Ensure predictions exist
  expect_true(all(c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w") %in% names(preds)))

  # Ensure predictions are numeric and non-NA
  for (col in c("y_hat_w", "y_hat_gc_w", "y_hat_gf_w")) {
    expect_false(any(is.na(preds[[col]])))
    expect_type(preds[[col]], "double")
  }
})

# ---- NUMERICAL ACCURACY TESTS (use large n) ----
test_that("hapr estimates coefficients correctly across parameter grid (probit)", {
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

    # Simulate dataset with current parameters
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
      model_type = "probit",
      improvement_ratio = 1 / (1 - params$var_epsilon)
    )

    # Extract estimated beta coefficients
    beta_hat <- fit$coefficients$beta
    err <- abs(beta_hat["w1"] - data$beta_w1)
    expect_lt(err, 0.07)

    expect_true("w2B" %in% names(beta_hat))
    expect_true("w2C" %in% names(beta_hat))
    expect_lt(abs(beta_hat["w2B"]), 0.07)
    expect_lt(abs(beta_hat["w2C"]), 0.07)
  }
})

# ---- SNAPSHOT TEST (use small n) ----
test_that("hapr print output is stable (probit)", {
  data <- simulate_mock_dataset(n = 100)

  # Fit the model
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "probit",
    improvement_ratio = 1.5
  )

  expect_snapshot(print(fit))
})

# ---- CI COVERAGE TEST ----
test_that("hapr confidence intervals for beta achieve exact match coverage (probit)", {
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

    latent <- beta_gf * gf + rnorm(n) + beta_w1 * w$w1
    y <- as.factor(as.numeric(latent > 0))

    list(
      y = y,
      gc = gc_normalized,
      w = w,
      gf = gf,
      true_improvement_ratio = true_improvement_ratio
    )
  }

  # ---- Part 2: True beta from simulation config ----
  beta_names <- c("gf", "w1", "w2B", "w2C")
  true_beta <- c(gf = default_params$beta_gf,
                 w1 = default_params$beta_w1,
                 w2B = 0,
                 w2C = 0)

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
      model_type = "probit",
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
    Coverage = round(colMeans(covered_matrix, na.rm = TRUE), 3)
  )

  print(coverage_df)

  expect_true(all(coverage_df$Coverage > 0.90),
              info = paste("Coverage below threshold:\n",
                           paste(coverage_df$Term, coverage_df$Coverage, collapse = "\n")))
})

# ---- PSEUDO-R2 CALCULATION TEST ----
test_that("pseudo-R2 from predict() matches manual calculation (probit)", {
  # Use existing helper function
  data <- simulate_mock_dataset(n = 1000)
  
  # Fit the hapr model
  fit <- hapr(
    y = data$y,
    gc = data$gc,
    w = data$w,
    model_type = "probit",
    r2_future = 0.31
  )

  # Simulate data
  sim_data <- hapr_simulate(fit, w = data$w, gc = data$gc)

  # Make predictions
  pred <- predict(fit, newdata = sim_data, covariates = "gf_w", type = "link")
  y_hat_gf_w <- pred$y_hat_gf_w

  # Calculate pseudo-R2 from predict function
  pseudo_r2_predict <- var(y_hat_gf_w) / (1 + var(y_hat_gf_w))

  # Manual calculation using the full formula
  beta_vec <- fit$coefficients$beta
  beta_w <- beta_vec[names(beta_vec) != "gf"]
  beta_gf <- beta_vec["gf"]

  # Get variance-covariance matrix of w
  w_matrix <- model.matrix(~., data = data$w)
  var_w <- var(w_matrix[, -1])  # Remove intercept

  # Calculate covariance between gf and w
  gf_sim <- sim_data$gf
  cov_gf_w <- cov(gf_sim, w_matrix[, -1])

  # Manual calculation
  beta_w_aligned <- beta_w[colnames(var_w)]
  cov_gf_w_col <- t(as.matrix(cov_gf_w))

  var_L_hat_gf_w <- 
    var(gf_sim) * beta_gf^2 +
    as.numeric(t(beta_w_aligned) %*% var_w %*% beta_w_aligned) + 
    2 * beta_gf * as.numeric(t(beta_w_aligned) %*% cov_gf_w_col)

  pseudo_r2_manual <- var_L_hat_gf_w / (1 + var_L_hat_gf_w)

  # Test that the two calculations match
  expect_equal(pseudo_r2_predict, as.numeric(pseudo_r2_manual), tolerance = 1e-10)
})
