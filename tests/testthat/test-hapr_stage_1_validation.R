test_that("hapr_first_stage rejects invalid model_type", {
  set.seed(123)
  
  # Create valid test data
  n <- 50
  w <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  gc <- rnorm(n)
  y <- rnorm(n)
  
  # Test invalid model_type
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w, model_type = "invalid"),
    "model_type must be one of: 'lm', 'probit'"
  )
  
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w, model_type = "cox"),
    "model_type must be one of: 'lm', 'probit'"
  )
})

test_that("hapr_first_stage rejects w as dataframe", {
  set.seed(123)
  
  n <- 50
  w_df <- data.frame(w1 = rnorm(n), w2 = rnorm(n))
  gc <- rnorm(n)
  y <- rnorm(n)
  
  # Should error because w must be a matrix
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w_df, model_type = "lm"),
    "w must be a numeric matrix"
  )
})

test_that("hapr_first_stage rejects w with constant column", {
  set.seed(123)
  
  n <- 50
  # Create w with a constant column
  w <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  w[, 2] <- 5  # Make second column constant
  gc <- rnorm(n)
  y <- rnorm(n)
  
  # Should error because w has a constant column
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w, model_type = "lm"),
    "w contains constant columns"
  )
  
  # Test with all zeros
  w_zero <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  w_zero[, 1] <- 0  # Make first column all zeros
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w_zero, model_type = "lm"),
    "w contains constant columns"
  )
})

test_that("hapr_first_stage rejects w that is not full rank", {
  set.seed(123)
  
  n <- 50
  # Create w with linearly dependent columns
  w1 <- rnorm(n)
  w <- cbind(w1, w1, rnorm(n))  # First two columns are identical
  gc <- rnorm(n)
  y <- rnorm(n)
  
  # Should error because w has linearly dependent columns
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w, model_type = "lm"),
    "w contains linearly dependent columns"
  )
  
  # Test with columns that are linear combinations
  w2 <- rnorm(n)
  w3 <- 2 * w1 + 3 * w2  # w3 is a linear combination of w1 and w2
  w_lin <- cbind(w1, w2, w3)
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w_lin, model_type = "lm"),
    "w contains linearly dependent columns"
  )
})

test_that("hapr_first_stage rejects empty w", {
  set.seed(123)
  
  n <- 50
  w_empty <- matrix(nrow = n, ncol = 0)
  gc <- rnorm(n)
  y <- rnorm(n)
  
  # Should error because w must have at least one column
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w_empty, model_type = "lm"),
    "w must have at least one column"
  )
})

test_that("hapr_first_stage rejects mismatched dimensions", {
  set.seed(123)
  
  n <- 50
  w <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  gc <- rnorm(n)
  y <- rnorm(n)
  
  # Test mismatched y length
  expect_error(
    hapr_first_stage(y = y[1:(n-1)], gc = gc, w = w, model_type = "lm"),
    "y, gc, and w must have the same number of observations"
  )
  
  # Test mismatched gc length
  expect_error(
    hapr_first_stage(y = y, gc = gc[1:(n-1)], w = w, model_type = "lm"),
    "y, gc, and w must have the same number of observations"
  )
  
  # Test mismatched w rows
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w[1:(n-1), , drop = FALSE], model_type = "lm"),
    "y, gc, and w must have the same number of observations"
  )
})

test_that("hapr_first_stage rejects missing values", {
  set.seed(123)
  
  n <- 50
  w <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  gc <- rnorm(n)
  y <- rnorm(n)
  
  # Test missing values in y
  y_na <- y
  y_na[1] <- NA
  expect_error(
    hapr_first_stage(y = y_na, gc = gc, w = w, model_type = "lm"),
    "y contains missing values"
  )
  
  # Test missing values in gc
  gc_na <- gc
  gc_na[1] <- NA
  expect_error(
    hapr_first_stage(y = y, gc = gc_na, w = w, model_type = "lm"),
    "gc contains missing values"
  )
  
  # Test missing values in w
  w_na <- w
  w_na[1, 1] <- NA
  expect_error(
    hapr_first_stage(y = y, gc = gc, w = w_na, model_type = "lm"),
    "w contains missing values"
  )
})

test_that("hapr_first_stage validates probit y is binary", {
  set.seed(123)
  
  n <- 50
  w <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  gc <- rnorm(n)
  
  # Test with non-binary y for probit
  y_continuous <- rnorm(n)
  expect_error(
    hapr_first_stage(y = y_continuous, gc = gc, w = w, model_type = "probit"),
    "For 'probit' model_type, y must be a binary numeric vector"
  )
  
  # Test with y containing values other than 0/1
  y_invalid <- c(rep(0, n-1), 2)
  expect_error(
    hapr_first_stage(y = y_invalid, gc = gc, w = w, model_type = "probit"),
    "For 'probit' model_type, y must be a binary numeric vector"
  )
})
