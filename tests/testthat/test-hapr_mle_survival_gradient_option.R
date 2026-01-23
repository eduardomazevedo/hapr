

test_that("Survival MLE runs without OpenMP", {
  set.seed(123)

  n <- 200
  var_epsilon <- 0.7
  var_v <- (1 - var_epsilon) * 0.4
  improvement_ratio <- 1 / (1 - var_epsilon)

  beta_g <- 0.6
  beta_w <- c(0.1, -0.2, 0.15, 0.05)
  theta <- c(0.0, 0.1, -0.25, 0.2)

  data <- mock_dataset_survival_exponential(
    n = n,
    var_v = var_v,
    var_epsilon = var_epsilon,
    beta_g = beta_g,
    beta_w = beta_w,
    theta = theta,
    censor_rate = 0.2
  )

  start_beta <- rep(0, ncol(data$w) + 2)

  fit <- hapr_mle_survival(
    event_time = data$event_time,
    event_status = data$event_status,
    gc = data$gc,
    w = data$w,
    improvement_ratio = improvement_ratio,
    model_type = "exponential",
    start_beta = start_beta,
    start_delta = numeric(0),
    use_openmp = FALSE,
    control = list(maxit = 100)
  )

  expect_s3_class(fit, "hapr_mle_fit")
  expect_true(all(is.finite(fit$parameters$beta)))
})
