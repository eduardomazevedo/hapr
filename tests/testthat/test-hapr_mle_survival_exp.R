#' Test point estimates for survival MLE across scenarios
#'
#' Tests that point estimates are within 4 standard errors of true coefficients
#' for various combinations of var_epsilon and n.
#' Produces artifact tables comparing true coefficients, estimates, SEs, and CIs.

test_that("Survival MLE point estimates are within 4 SE", {
  mock_dataset_path <- file.path("dev", "mock_dataset.R")
  if (file.exists(mock_dataset_path)) {
    source(mock_dataset_path)
  } else {
    mock_dataset_path <- file.path("..", "..", "dev", "mock_dataset.R")
    if (file.exists(mock_dataset_path)) {
      source(mock_dataset_path)
    } else {
      stop("Cannot find dev/mock_dataset.R")
    }
  }

  BETA_G <- 0.6
  BETA_W <- c(0.1, -0.2, 0.15, 0.05)
  THETA <- c(0.0, 0.1, -0.25, 0.2)
  CENSOR_RATE <- 0.2

  VAR_EPSILON_VALUES <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  N_VALUES <- c(1e3, 1e4)
  MODEL_TYPES <- c("exponential", "weibull")
  LOG_K_VALUES <- c(-1, 0, 1)

  artifact_dir <- testthat::test_path("_artifacts", "mle_survival")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }

  all_results <- list()

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
          runtime_mle_ms <- mle_time[["elapsed"]] * 1000

          true_coef <- c(BETA_G, BETA_W)
          names(true_coef) <- c("gf", "(Intercept)", "w1", "w2", "w3")

          estimates <- mle_fit$parameters$beta[names(true_coef)]
          se <- if (is.null(mle_fit$standard_errors)) {
            rep(NA_real_, length(true_coef))
          } else {
            mle_fit$standard_errors[names(true_coef)]
          }
          names(se) <- names(true_coef)

          lower_ci <- estimates - 1.96 * se
          upper_ci <- estimates + 1.96 * se

          differences <- abs(estimates - true_coef)
          within_4se <- differences <= 4 * se

          comparison_table <- data.frame(
            Coefficient = names(true_coef),
            True_Value = true_coef,
            Point_Estimate = estimates,
            Std_Error = se,
            Lower_CI = lower_ci,
            Upper_CI = upper_ci,
            Difference = differences,
            Within_4SE = within_4se,
            Runtime_MLE_ms = rep(runtime_mle_ms, length(true_coef)),
            row.names = NULL,
            stringsAsFactors = FALSE
          )

          all_results[[scenario_name]] <- list(
            scenario = scenario_name,
            var_epsilon = var_epsilon,
            n = n,
            model_type = model_type,
            log_k = log_k,
            comparison_table = comparison_table,
            all_within_4se = all(within_4se, na.rm = TRUE)
          )

          expect_true(
            all(within_4se, na.rm = TRUE),
            info = sprintf(
              "Scenario %s: Some estimates are not within 4 SE of true values.\n%s",
              scenario_name,
              paste(capture.output(print(comparison_table[!within_4se & !is.na(within_4se), ])),
                    collapse = "\n")
            )
          )

          artifact_file <- file.path(artifact_dir, paste0(scenario_name, ".csv"))
          write.csv(comparison_table, artifact_file, row.names = FALSE)
        }
      }
    }
  }

  summary_table <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      n = x$n,
      var_epsilon = x$var_epsilon,
      model_type = x$model_type,
      log_k = x$log_k,
      All_Within_4SE = x$all_within_4se,
      stringsAsFactors = FALSE
    )
  }))

  summary_file <- file.path(artifact_dir, "summary.csv")
  write.csv(summary_table, summary_file, row.names = FALSE)

  expect_true(
    all(sapply(all_results, function(x) x$all_within_4se), na.rm = TRUE),
    info = sprintf(
      "Not all scenarios passed. Summary:\n%s",
      paste(capture.output(print(summary_table)), collapse = "\n")
    )
  )
})
