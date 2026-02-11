#' Test point estimates for survival MLE across scenarios
#'
#' Tests that point estimates are within 4 standard errors of true coefficients
#' for various combinations of var_epsilon and n.
#' Produces artifact tables comparing true coefficients, estimates, SEs, and CIs.

test_that("Survival MLE point estimates are within 4 SE", {
  params <- survival_params_default()
  scenarios <- survival_scenarios(is_slow_enabled(), log_k_values = c(-1, 0, 1))
  artifact_dir <- ensure_artifact_dir("mle_survival")

  all_results <- run_survival_tests(
    test_type = "point",
    params = params,
    scenarios = scenarios,
    artifact_dir = artifact_dir
  )

  summary_table <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      n = x$n,
      var_epsilon = x$var_epsilon,
      model_type = x$model_type,
      log_k = x$log_k,
      var_v_factor = x$var_v_factor,
      All_Within_4SE = x$all_within_4se,
      stringsAsFactors = FALSE
    )
  }))

  write_summary_csv(artifact_dir, "summary.csv", summary_table)

  expect_true(
    all(sapply(all_results, function(x) x$all_within_4se), na.rm = TRUE),
    info = sprintf(
      "Not all scenarios passed. Summary:\n%s",
      paste(capture.output(print(summary_table)), collapse = "\n")
    )
  )
})


test_that("Survival MLE runs without OpenMP", {
  set.seed(123)

  n <- 200
  var_epsilon <- 0.7
  var_v <- (1 - var_epsilon) * 0.4
  improvement_ratio <- 1 / (1 - var_epsilon)

  beta_g <- 0.6
  beta_w <- c(0.1, -0.2, 0.15, 0.05)
  theta <- normalize_theta(
    theta = c(0.0, 0.1, -0.25, 0.2),
    var_v = var_v,
    var_epsilon = var_epsilon
  )

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
