#' Test survival MLE estimator with analytic gradients
#'
#' Compares estimates against the existing MLE estimator and true parameters.
#' Produces artifacts with parameter diffs and runtimes.

test_that("Survival MLE gradient estimator matches existing estimator", {
  BETA_G <- 0.6
  BETA_W <- c(0.1, -0.2, 0.15, 0.05)
  THETA <- c(0.0, 0.1, -0.25, 0.2)
  CENSOR_RATE <- 0.2

  VAR_EPSILON_VALUES <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  run_slow_env <- Sys.getenv("RUN_SLOW_TESTS", unset = "false")
  run_slow <- tolower(run_slow_env) %in% c("true", "1", "yes")
  N_VALUES <- if (run_slow) c(1e3, 1e4, 1e5) else c(1e3, 1e4)
  MODEL_TYPES <- c("exponential", "weibull")
  LOG_K_VALUES <- c(-1, 0, 1)

  artifact_dir <- testthat::test_path("_artifacts", "mle_survival_grad")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }

  results <- list()
  summary_results <- list()

  for (var_epsilon in VAR_EPSILON_VALUES) {
    for (n in N_VALUES) {
      for (model_type in MODEL_TYPES) {
        log_k_values <- if (model_type == "weibull") LOG_K_VALUES else 0
        for (log_k in log_k_values) {
          scenario_name <- sprintf("%s_n%d_ve%.1f_lk%.1f",
                                   substr(model_type, 1, 3),
                                   as.integer(n),
                                   var_epsilon,
                                   log_k)

          set.seed(123 + which(VAR_EPSILON_VALUES == var_epsilon) * 100 +
                     which(N_VALUES == n) * 10 +
                     which(MODEL_TYPES == model_type) * 1000 +
                     round(log_k * 10))

          var_v <- (1 - var_epsilon) * 0.4
          improvement_ratio <- 1 / (1 - var_epsilon)

          if (model_type == "exponential") {
            data <- mock_dataset_survival_exponential(
              n = n,
              var_v = var_v,
              var_epsilon = var_epsilon,
              beta_g = BETA_G,
              beta_w = BETA_W,
              theta = THETA,
              censor_rate = CENSOR_RATE
            )
          } else {
            data <- mock_dataset_survival_weibull(
              n = n,
              var_v = var_v,
              var_epsilon = var_epsilon,
              beta_g = BETA_G,
              beta_w = BETA_W,
              theta = THETA,
              log_k = log_k,
              censor_rate = CENSOR_RATE
            )
          }

          start_beta <- rep(0, ncol(data$w) + 2)
          start_delta <- if (model_type == "weibull") c(log_k = 0) else numeric(0)

          mle_time <- system.time({
            mle_fit <- hapr_mle_survival(
              event_time = data$event_time,
              event_status = data$event_status,
              gc = data$gc,
              w = data$w,
              improvement_ratio = improvement_ratio,
              model_type = model_type,
              start_beta = start_beta,
              start_delta = start_delta,
              control = list(maxit = 150)
            )
          })

          grad_time <- system.time({
            grad_fit <- hapr_mle_survival_grad(
              event_time = data$event_time,
              event_status = data$event_status,
              gc = data$gc,
              w = data$w,
              improvement_ratio = improvement_ratio,
              model_type = model_type,
              start_beta = start_beta,
              start_delta = start_delta,
              control = list(maxit = 150)
            )
          })

          runtime_mle_ms <- mle_time[["elapsed"]] * 1000
          runtime_grad_ms <- grad_time[["elapsed"]] * 1000

          true_coef <- c(BETA_G, BETA_W)
          names(true_coef) <- c("gf", "(Intercept)", "w1", "w2", "w3")

          est_mle <- mle_fit$parameters$beta[names(true_coef)]
          est_grad <- grad_fit$parameters$beta[names(true_coef)]

          diff_grad_vs_mle <- est_grad - est_mle
          diff_grad_vs_true <- est_grad - true_coef

          within_4se <- rep(NA_real_, length(true_coef))
          if (!is.null(grad_fit$standard_errors)) {
            se <- grad_fit$standard_errors[names(true_coef)]
            within_4se <- abs(diff_grad_vs_true) <= 4 * se
          }

          results[[scenario_name]] <- data.frame(
            Scenario = scenario_name,
            Coefficient = names(true_coef),
            Estimate_MLE = est_mle,
            Estimate_Grad = est_grad,
            Diff_Grad_vs_MLE = diff_grad_vs_mle,
            Diff_Grad_vs_True = diff_grad_vs_true,
            Within_4SE = within_4se,
            Runtime_MLE_ms = rep(runtime_mle_ms, length(true_coef)),
            Runtime_Grad_ms = rep(runtime_grad_ms, length(true_coef)),
            stringsAsFactors = FALSE
          )

          summary_results[[scenario_name]] <- data.frame(
            Scenario = scenario_name,
            n = n,
            Var_Epsilon = var_epsilon,
            Model = model_type,
            Log_k = log_k,
            Max_Abs_Diff_Grad_vs_MLE = max(abs(diff_grad_vs_mle)),
            All_Within_4SE = all(within_4se, na.rm = TRUE),
            Runtime_MLE_ms = runtime_mle_ms,
            Runtime_Grad_ms = runtime_grad_ms,
            stringsAsFactors = FALSE
          )

          expect_true(max(abs(diff_grad_vs_mle)) < 1e-3)
          if (all(!is.na(within_4se))) {
            expect_true(all(within_4se))
          }
        }
      }
    }
  }

  detail_table <- do.call(rbind, results)
  summary_table <- do.call(rbind, summary_results)
  write.csv(detail_table,
            file.path(artifact_dir, "estimator_comparison_detail.csv"),
            row.names = FALSE)
  write.csv(summary_table,
            file.path(artifact_dir, "estimator_comparison_summary.csv"),
            row.names = FALSE)
})
