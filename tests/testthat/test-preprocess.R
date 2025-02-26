library(testthat)
library(survival)

test_that("preprocess works for linear models", {
  set.seed(123)
  n <- 100
  y <- rnorm(n)
  gc <- rnorm(n)
  w <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- preprocess(y, gc, w, model_type = "lm")

  expect_type(result$y, "double")
  expect_equal(length(result$y), n)
  expect_equal(length(result$gc), n)
  expect_equal(nrow(result$w), n)
})

test_that("preprocess works for probit models with numeric y", {
  set.seed(123)
  n <- 100
  y <- rbinom(n, 1, 0.5)
  gc <- rnorm(n)
  w <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- preprocess(y, gc, w, model_type = "probit")

  expect_true(is.factor(result$y))
  expect_equal(nlevels(result$y), 2)
})

test_that("preprocess works for probit models with logical y", {
  set.seed(123)
  n <- 100
  y <- sample(c(TRUE, FALSE), n, replace = TRUE)
  gc <- rnorm(n)
  w <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- preprocess(y, gc, w, model_type = "probit")

  expect_true(is.factor(result$y))
  expect_equal(nlevels(result$y), 2)
})

test_that("preprocess works for probit models with factor y", {
  set.seed(123)
  n <- 100
  y <- factor(sample(c("A", "B"), n, replace = TRUE))
  gc <- rnorm(n)
  w <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- preprocess(y, gc, w, model_type = "probit")

  expect_true(is.factor(result$y))
  expect_equal(nlevels(result$y), 2)
})

test_that("preprocess works for cox models", {
  set.seed(123)
  n <- 100
  time <- rexp(n)
  status <- sample(0:1, n, replace = TRUE)
  y <- Surv(time, status)
  gc <- rnorm(n)
  w <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- preprocess(y, gc, w, model_type = "cox")

  expect_true(inherits(result$y, "Surv"))
  expect_equal(length(result$y), n)
})

test_that("preprocess normalizes gc", {
  set.seed(123)
  n <- 100
  y <- rnorm(n)
  gc <- rnorm(n)
  w <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

  result <- preprocess(y, gc, w, model_type = "lm")

  expect_true(abs(mean(result$gc)) < 1e-6)
  expect_true(abs(sd(result$gc) - 1) < 1e-6)
})

test_that("preprocess handles missing values", {
  set.seed(123)
  n <- 100
  y <- rnorm(n)
  gc <- rnorm(n)
  w <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  
  y[1:5] <- NA
  gc[6:10] <- NA
  w$x1[11:15] <- NA

  result <- preprocess(y, gc, w, model_type = "lm")

  expect_equal(length(result$y), n - 15)
  expect_equal(length(result$gc), n - 15)
  expect_equal(nrow(result$w), n - 15)
})

test_that("preprocess throws an error for invalid model_type", {
  set.seed(123)
  y <- rnorm(100)
  gc <- rnorm(100)
  w <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

  expect_error(preprocess(y, gc, w, model_type = "invalid"),
               "model_type must be one of: 'lm', 'probit', 'cox'")
})

test_that("preprocess throws an error for invalid probit y", {
  set.seed(123)
  y <- factor(c("A", "B", "C"))  # More than 2 levels
  gc <- rnorm(100)
  w <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

  expect_error(preprocess(y, gc, w, model_type = "probit"),
               "For probit models, y must have exactly 2 levels")
})

test_that("preprocess throws an error for non-Surv cox model y", {
  set.seed(123)
  y <- rnorm(100)  # Not a Surv object
  gc <- rnorm(100)
  w <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

  expect_error(preprocess(y, gc, w, model_type = "cox"),
               "For cox models, y must be a Surv object")
})

test_that("preprocess throws an error when lengths do not match", {
  set.seed(123)
  y <- rnorm(100)
  gc <- rnorm(101)  # Different length
  w <- data.frame(x1 = rnorm(100), x2 = rnorm(100))

  expect_error(preprocess(y, gc, w, model_type = "lm"),
               "y, gc, and w must have the same number of observations")
})

test_that("preprocess throws an error if w is not a data frame", {
  set.seed(123)
  y <- rnorm(100)
  gc <- rnorm(100)
  w <- matrix(rnorm(200), ncol = 2)  # Not a data frame

  expect_error(preprocess(y, gc, w, model_type = "lm"),
               "w must be a data frame")
})
