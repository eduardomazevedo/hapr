test_that("fit_lm matches lm output for linear regression with intercept", {
  set.seed(123)
  
  # Simulate data
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- c("x1", "x2", "x3")
  beta_true <- c(2, -1, 0.5, 1.5)  # intercept + 3 predictors
  y <- beta_true[1] + X %*% beta_true[-1] + rnorm(n, sd = 0.5)
  
  # Fit with low-level function
  fit_low <- hapr:::fit_lm(y = y, X = X, add_intercept = TRUE)
  
  # Fit with standard lm
  df_lm <- data.frame(y = y, X)
  fit_std <- lm(y ~ ., data = df_lm)
  
  # Compare coefficients (should match exactly)
  expect_equal(
    fit_low$coefficients,
    coef(fit_std),
    tolerance = 1e-10
  )
  
  # Compare vcov (should match exactly)
  expect_equal(
    fit_low$vcov_coefficients,
    vcov(fit_std),
    tolerance = 1e-10
  )
  
  # Compare sigma squared
  expect_equal(
    fit_low$sigma_squared,
    summary(fit_std)$sigma^2,
    tolerance = 1e-10
  )
  
  # Compare R-squared
  expect_equal(
    fit_low$r2,
    summary(fit_std)$r.squared,
    tolerance = 1e-10
  )
})

test_that("fit_lm matches lm output for linear regression without intercept", {
  set.seed(123)
  
  # Simulate data
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- c("x1", "x2", "x3")
  beta_true <- c(2, -1, 0.5, 1.5)
  y <- beta_true[1] + X %*% beta_true[-1] + rnorm(n, sd = 0.5)
  
  # Fit without intercept
  fit_low <- hapr:::fit_lm(y = y, X = X, add_intercept = FALSE)
  df_lm <- data.frame(y = y, X)
  fit_std <- lm(y ~ 0 + ., data = df_lm)
  
  # Compare coefficients
  expect_equal(
    fit_low$coefficients,
    coef(fit_std),
    tolerance = 1e-10
  )
  
  # Compare vcov
  expect_equal(
    fit_low$vcov_coefficients,
    vcov(fit_std),
    tolerance = 1e-10
  )
  
  # Compare sigma squared
  expect_equal(
    fit_low$sigma_squared,
    summary(fit_std)$sigma^2,
    tolerance = 1e-10
  )
})

test_that("fit_probit matches glm output for probit regression", {
  set.seed(123)
  
  # Simulate data for probit regression
  n <- 200
  p <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- c("x1", "x2")
  beta_true <- c(-0.5, 1, -0.8)  # intercept + 2 predictors
  linear_pred <- beta_true[1] + X %*% beta_true[-1]
  prob <- pnorm(linear_pred)
  y <- rbinom(n, size = 1, prob = prob)
  
  # Fit with low-level function
  fit_low <- hapr:::fit_probit(y = y, X = X, add_intercept = TRUE)
  
  # Fit with standard glm
  df_probit <- data.frame(y = y, X)
  fit_std <- glm(y ~ ., data = df_probit, family = binomial(link = "probit"))
  
  # Compare coefficients (should match exactly)
  expect_equal(
    fit_low$coefficients,
    coef(fit_std),
    tolerance = 1e-10
  )
  
  # Compare vcov (should match exactly)
  expect_equal(
    fit_low$vcov_coefficients,
    vcov(fit_std),
    tolerance = 1e-10
  )
  
  # R-squared uses different definitions, so we just check it's a valid value
  expect_true(fit_low$r2 >= 0 && fit_low$r2 <= 1)
})

test_that("fit_probit matches glm output for probit regression without intercept", {
  set.seed(123)
  
  # Simulate data for probit regression
  n <- 200
  p <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- c("x1", "x2")
  beta_true <- c(1, -0.8)  # no intercept, 2 predictors
  linear_pred <- X %*% beta_true
  prob <- pnorm(linear_pred)
  y <- rbinom(n, size = 1, prob = prob)
  
  # Fit without intercept
  fit_low <- hapr:::fit_probit(y = y, X = X, add_intercept = FALSE)
  df_probit <- data.frame(y = y, X)
  fit_std <- glm(y ~ 0 + ., data = df_probit, family = binomial(link = "probit"))
  
  # Compare coefficients
  expect_equal(
    fit_low$coefficients,
    coef(fit_std),
    tolerance = 1e-10
  )
  
  # Compare vcov
  expect_equal(
    fit_low$vcov_coefficients,
    vcov(fit_std),
    tolerance = 1e-10
  )
})

test_that("fit_lm returns correct structure and types", {
  set.seed(123)
  
  n <- 50
  p <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- rnorm(n)
  
  result <- hapr:::fit_lm(y = y, X = X, add_intercept = TRUE)
  
  # Check structure
  expect_named(result, c("coefficients", "vcov_coefficients", "sigma_squared", 
                         "var_sigma_squared", "r2"))
  
  # Check types
  expect_type(result$coefficients, "double")
  expect_type(result$vcov_coefficients, "double")
  expect_type(result$sigma_squared, "double")
  expect_type(result$var_sigma_squared, "double")
  expect_type(result$r2, "double")
  
  # Check dimensions
  expect_length(result$coefficients, p + 1)  # intercept + p predictors
  expect_equal(dim(result$vcov_coefficients), c(p + 1, p + 1))
  
  # Check valid ranges
  expect_true(result$sigma_squared > 0)
  expect_true(result$var_sigma_squared > 0)
  expect_true(result$r2 >= 0 && result$r2 <= 1)
})

test_that("fit_probit returns correct structure and types", {
  set.seed(123)
  
  n <- 100
  p <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- rbinom(n, size = 1, prob = 0.5)
  
  result <- hapr:::fit_probit(y = y, X = X, add_intercept = TRUE)
  
  # Check structure
  expect_named(result, c("coefficients", "vcov_coefficients", "r2"))
  
  # Check types
  expect_type(result$coefficients, "double")
  expect_type(result$vcov_coefficients, "double")
  expect_type(result$r2, "double")
  
  # Check dimensions
  expect_length(result$coefficients, p + 1)  # intercept + p predictors
  expect_equal(dim(result$vcov_coefficients), c(p + 1, p + 1))
  
  # Check valid ranges
  expect_true(result$r2 >= 0 && result$r2 <= 1)
})

test_that("fit_lm handles edge case with single predictor", {
  set.seed(123)
  
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "x1"
  y <- 2 + 1.5 * X[, 1] + rnorm(n, sd = 0.5)
  
  fit_low <- hapr:::fit_lm(y = y, X = X, add_intercept = TRUE)
  df_lm <- data.frame(y = y, x1 = X[, 1])
  fit_std <- lm(y ~ ., data = df_lm)
  
  expect_equal(fit_low$coefficients, coef(fit_std), tolerance = 1e-10)
  expect_equal(fit_low$vcov_coefficients, vcov(fit_std), tolerance = 1e-10)
})

test_that("fit_probit handles edge case with single predictor", {
  set.seed(123)
  
  n <- 100
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "x1"
  linear_pred <- -0.5 + 1 * X[, 1]
  prob <- pnorm(linear_pred)
  y <- rbinom(n, size = 1, prob = prob)
  
  fit_low <- hapr:::fit_probit(y = y, X = X, add_intercept = TRUE)
  df_probit <- data.frame(y = y, x1 = X[, 1])
  fit_std <- glm(y ~ ., data = df_probit, family = binomial(link = "probit"))
  
  expect_equal(fit_low$coefficients, coef(fit_std), tolerance = 1e-10)
  expect_equal(fit_low$vcov_coefficients, vcov(fit_std), tolerance = 1e-10)
})
