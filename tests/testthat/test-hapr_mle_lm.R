#' Test MLE estimator against two-stage linear model

test_that("hapr_mle matches two-stage estimates for linear model", {
  source(testthat::test_path("..", "..", "dev", "mock_dataset.R"))

  var_epsilon_values <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  n_values <- c(1e3, 1e4, 1e5)

  beta_g <- 1.2
  beta_w <- c(0.2, -0.1, 0.3, -0.2)
  theta <- c(0.0, 0.1, -0.2, 0.25)
  var_y <- 1.0

  loglik_fn <- function(y, linpred, delta) {
    sigma <- exp(delta["log_sigma"])
    stats::dnorm(y, mean = linpred, sd = sigma, log = TRUE)
  }

  artifact_dir <- testthat::test_path("_artifacts", "mle_lm")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }

  for (var_epsilon in var_epsilon_values) {
    for (n in n_values) {
      set.seed(123 + which(var_epsilon_values == var_epsilon) * 100 +
                 which(n_values == n) * 10)

      var_v <- (1 - var_epsilon) * 0.5
      improvement_ratio <- 1 / (1 - var_epsilon)

      data <- mock_dataset_lm(
        n = n,
        var_v = var_v,
        var_epsilon = var_epsilon,
        beta_g = beta_g,
        beta_w = beta_w,
        theta = theta,
        var_y = var_y
      )

      first_stage <- hapr_first_stage(
        y = data$y,
        gc = data$gc,
        w = data$w,
        model_type = "lm"
      )
      two_stage_time <- system.time({
        first_stage <- hapr_first_stage(
          y = data$y,
          gc = data$gc,
          w = data$w,
          model_type = "lm"
        )
        second_stage <- hapr_second_stage(
          first_stage = first_stage,
          improvement_ratio = improvement_ratio
        )
      })

      n_params <- ncol(data$w) + 2
      start_beta <- rep(0, n_params)
      start_delta <- c(log_sigma = 0.0)

      mle_time <- system.time({
        mle_fit <- hapr_mle(
          y = data$y,
          gc = data$gc,
          w = data$w,
          improvement_ratio = improvement_ratio,
          loglik_fn = loglik_fn,
          start_beta = start_beta,
          start_delta = start_delta,
          control = list(maxit = 150)
        )
      })
      runtime_two_stage_ms <- two_stage_time[["elapsed"]] * 1000
      runtime_mle_ms <- mle_time[["elapsed"]] * 1000

      true_coef <- c(beta_g, beta_w)
      names(true_coef) <- c("gf", "(Intercept)", "w1", "w2", "w3")

      mle_beta <- mle_fit$parameters$beta[names(true_coef)]
      two_stage_beta <- second_stage$parameters$beta[names(true_coef)]
      two_stage_se <- second_stage$standard_errors[names(true_coef)]
      if (is.null(mle_fit$standard_errors)) {
        mle_se <- rep(NA_real_, length(true_coef))
        names(mle_se) <- names(true_coef)
      } else {
        mle_se <- mle_fit$standard_errors[names(true_coef)]
      }

      comparison_table <- data.frame(
        Parameter = names(true_coef),
        True_Value = true_coef,
        MLE_Estimate = mle_beta,
        Two_Stage_Estimate = two_stage_beta,
        Two_Stage_SE = two_stage_se,
        MLE_SE = mle_se,
        SE_Ratio = mle_se / two_stage_se,
        Runtime_Two_Stage_ms = rep(runtime_two_stage_ms, length(true_coef)),
        Runtime_MLE_ms = rep(runtime_mle_ms, length(true_coef)),
        row.names = NULL,
        stringsAsFactors = FALSE
      )

      scenario_name <- sprintf("lm_n%d_ve%.1f", as.integer(n), var_epsilon)
      artifact_file <- file.path(artifact_dir, paste0(scenario_name, ".csv"))
      write.csv(comparison_table, artifact_file, row.names = FALSE)

      expect_equal(
        mle_fit$parameters$beta,
        second_stage$parameters$beta,
        tolerance = 5e-2,
        info = sprintf("var_epsilon=%.1f n=%d", var_epsilon, n)
      )
      expect_equal(mle_fit$opt$convergence, 0)
    }
  }
})
