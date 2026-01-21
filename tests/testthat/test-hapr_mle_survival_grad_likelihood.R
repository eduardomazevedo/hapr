#' Test analytic gradients for survival MLE likelihood
#'
#' Compares analytic gradients against finite-difference gradients and
#' verifies NLL matches the existing likelihood implementation.
#' Produces artifacts with gradient errors and timings.

test_that("Survival MLE analytic gradients match finite differences", {
  set.seed(123)

  BETA_G <- 0.4
  BETA_W <- c(0.05, -0.1, 0.12, 0.03)
  THETA <- c(0.0, 0.08, -0.12, 0.1)
  LOG_K <- 0.3
  CENSOR_RATE <- 0.2

  scenarios <- list(
    list(model_type = "exponential", n = 300, var_epsilon = 0.7, log_k = 0),
    list(model_type = "weibull", n = 300, var_epsilon = 0.7, log_k = LOG_K)
  )

  artifact_dir <- testthat::test_path("_artifacts", "mle_survival_grad")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }

  finite_diff_grad <- function(fn, params, eps = 1e-6) {
    grad <- numeric(length(params))
    for (i in seq_along(params)) {
      delta <- rep(0, length(params))
      delta[i] <- eps
      grad[i] <- (fn(params + delta) - fn(params - delta)) / (2 * eps)
    }
    grad
  }

  results <- list()

  for (scenario in scenarios) {
    var_epsilon <- scenario$var_epsilon
    var_v <- (1 - var_epsilon) * 0.4

    if (scenario$model_type == "exponential") {
      data <- mock_dataset_survival_exponential(
        n = scenario$n,
        var_v = var_v,
        var_epsilon = var_epsilon,
        beta_g = BETA_G,
        beta_w = BETA_W,
        theta = THETA,
        censor_rate = CENSOR_RATE
      )
    } else {
      data <- mock_dataset_survival_weibull(
        n = scenario$n,
        var_v = var_v,
        var_epsilon = var_epsilon,
        beta_g = BETA_G,
        beta_w = BETA_W,
        theta = THETA,
        log_k = scenario$log_k,
        censor_rate = CENSOR_RATE
      )
    }

    first_stage <- hapr_first_stage(
      y = data$event_time,
      gc = data$gc,
      w = data$w,
      model_type = "mle"
    )

    gc <- first_stage$preprocessed$gc
    w <- first_stage$preprocessed$w
    X_w <- cbind(1, w)
    colnames(X_w)[1] <- "(Intercept)"
    w_theta <- c(X_w %*% first_stage$parameters$theta)
    posterior <- abc(var_epsilon, var_v)

    nll_old <- make_hapr_mle_likelihood_survival(
      event_time = data$event_time,
      event_status = data$event_status,
      gc = gc,
      w_theta = w_theta,
      X_w = X_w,
      posterior = posterior,
      model_type = scenario$model_type
    )
    nll_grad <- make_hapr_mle_likelihood_survival_grad(
      event_time = data$event_time,
      event_status = data$event_status,
      gc = gc,
      w_theta = w_theta,
      X_w = X_w,
      posterior = posterior,
      model_type = scenario$model_type
    )

    params <- c(BETA_G, BETA_W)
    if (scenario$model_type == "weibull") {
      params <- c(params, log_k = scenario$log_k)
    }

    nll_old_val <- nll_old(params)
    nll_new_val <- nll_grad$fn(params)
    expect_equal(nll_new_val, nll_old_val, tolerance = 1e-8)

    analytic_grad <- nll_grad$gr(params)
    numeric_grad <- finite_diff_grad(nll_grad$fn, params)

    grad_diff <- analytic_grad - numeric_grad
    max_abs_diff <- max(abs(grad_diff))
    mean_abs_diff <- mean(abs(grad_diff))
    nll_diff <- abs(nll_new_val - nll_old_val)

    evals <- 20
    runtime_nll_ms <- system.time({
      for (i in seq_len(evals)) {
        nll_grad$fn(params)
      }
    })[["elapsed"]] * 1000 / evals

    runtime_grad_ms <- system.time({
      for (i in seq_len(evals)) {
        nll_grad$gr(params)
      }
    })[["elapsed"]] * 1000 / evals

    scenario_name <- sprintf("%s_n%d_ve%.1f",
                             substr(scenario$model_type, 1, 3),
                             scenario$n,
                             var_epsilon)

    results[[scenario_name]] <- data.frame(
      Scenario = scenario_name,
      Model = scenario$model_type,
      n = scenario$n,
      Var_Epsilon = var_epsilon,
      Max_Abs_Grad_Diff = max_abs_diff,
      Mean_Abs_Grad_Diff = mean_abs_diff,
      NLL_Diff = nll_diff,
      Runtime_NLL_ms = runtime_nll_ms,
      Runtime_Grad_ms = runtime_grad_ms,
      stringsAsFactors = FALSE
    )

    expect_true(max_abs_diff < 1e-4)
  }

  summary_table <- do.call(rbind, results)
  artifact_file <- file.path(artifact_dir, "likelihood_gradient_summary.csv")
  write.csv(summary_table, artifact_file, row.names = FALSE)
})
