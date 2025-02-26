library(testthat)

test_that("r2_liability_probit works correctly for a valid probit model", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  y <- rbinom(n, 1, pnorm(0.5 * x))  # Generate binary response using probit link
  
  model <- glm(y ~ x, family = binomial(link = "probit"))
  r2 <- r2_liability_probit(model)
  
  expect_type(r2, "double")
  expect_true(r2 >= 0 && r2 <= 1)
})

test_that("r2_liability_probit returns 0 for intercept-only model", {
  set.seed(123)
  n <- 100
  y <- rbinom(n, 1, 0.5)  # Random binary response with no predictors
  
  model <- glm(y ~ 1, family = binomial(link = "probit"))
  r2 <- r2_liability_probit(model)
  
  expect_equal(r2, 0)
})

test_that("r2_liability_probit throws an error for non-glm objects", {
  model <- lm(mpg ~ hp, data = mtcars)  # Not a glm object
  
  expect_error(r2_liability_probit(model), "Error: The input model is not a glm object.")
})

test_that("r2_liability_probit throws an error for non-probit link functions", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.5 * x))  # Logistic link instead of probit
  
  model <- glm(y ~ x, family = binomial(link = "logit"))
  
  expect_error(r2_liability_probit(model), "Error: The model does not use a probit link function.")
})