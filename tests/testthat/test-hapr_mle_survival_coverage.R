#' Coverage tests for survival MLE
#'
#' Runs multiple simulations per scenario and checks that coverage is above 85%
#' for all beta parameters. This test is slow and only runs when
#' RUN_SLOW_TESTS is set.

test_that("Survival MLE coverage intervals are above 85%", {
  run_slow <- Sys.getenv("RUN_SLOW_TESTS", unset = "false")
  if (!tolower(run_slow) %in% c("true", "1", "yes")) {
    skip("Coverage tests skipped. Set RUN_SLOW_TESTS=true to run.")
  }

  BETA_G <- 0.6
  BETA_W <- c(0.1, -0.2, 0.15, 0.05)
  THETA <- c(0.0, 0.1, -0.25, 0.2)
  CENSOR_RATE <- 0.2

  VAR_EPSILON_VALUES <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  N_VALUES <- c(1e3, 1e4, 1e5)
  MODEL_TYPES <- c("exponential", "weibull")
  LOG_K_VALUES <- c(0)
  N_SIMULATIONS <- 100

  artifact_dir <- testthat::test_path("_artifacts", "coverage_survival")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }

  all_results <- list()

  for (var_epsilon in VAR_EPSILON_VALUES) {
    for (n in N_VALUES) {
      for (model_type in MODEL_TYPES) {
        log_k_values <- if (model_type == "weibull") LOG_K_VALUES else 0
        for (log_k in log_k_values) {
          improvement_ratio <- 1 / (1 - var_epsilon)
          var_v <- (1 - var_epsilon) * 0.4

          scenario_name <- sprintf("%s_n%d_ve%.1f_lk%.1f",
                                   substr(model_type, 1, 3),
                                   as.integer(n),
                                   var_epsilon,
                                   log_k)

          true_coef <- c(BETA_G, BETA_W)
          names(true_coef) <- c("gf", "(Intercept)", "w1", "w2", "w3")
          n_coef <- length(true_coef)

          all_estimates <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
          colnames(all_estimates) <- names(true_coef)
          all_se <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
          colnames(all_se) <- names(true_coef)
          all_lower_ci <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
          colnames(all_lower_ci) <- names(true_coef)
          all_upper_ci <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
          colnames(all_upper_ci) <- names(true_coef)

          adjusted_count <- 0
          used_count <- 0
          for (sim in 1:N_SIMULATIONS) {
            set.seed(12345 +
                     which(VAR_EPSILON_VALUES == var_epsilon) * 10000 +
                     which(N_VALUES == n) * 1000 +
                     which(MODEL_TYPES == model_type) * 100 +
                     round(log_k * 10) +
                     sim)

            data <- if (model_type == "exponential") {
              mock_dataset_survival_exponential(
                n = n,
                var_v = var_v,
                var_epsilon = var_epsilon,
                beta_g = BETA_G,
                beta_w = BETA_W,
                theta = THETA,
                censor_rate = CENSOR_RATE
              )
            } else {
              mock_dataset_survival_weibull(
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

            first_stage <- hapr_first_stage(
              y = data$event_time,
              gc = data$gc,
              w = data$w,
              model_type = "mle"
            )
            max_improvement_ratio <- first_stage$stats$max_improvement_ratio
            if (n == 1e3 && improvement_ratio >= max_improvement_ratio) {
              adjusted_count <- adjusted_count + 1
              next
            }

            fit <- hapr_mle_survival(
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

            ci_beta <- fit$ci_beta
            if (!is.null(ci_beta)) {
              coef_names <- rownames(ci_beta)
              estimates <- setNames(ci_beta$Estimate, coef_names)
              se <- setNames(ci_beta$Std.Error, coef_names)
              lower_ci <- setNames(ci_beta$Lower, coef_names)
              upper_ci <- setNames(ci_beta$Upper, coef_names)

              all_estimates[sim, ] <- estimates[names(true_coef)]
              all_se[sim, ] <- se[names(true_coef)]
              all_lower_ci[sim, ] <- lower_ci[names(true_coef)]
              all_upper_ci[sim, ] <- upper_ci[names(true_coef)]
            }
            used_count <- used_count + 1
          }

          mean_estimates <- colMeans(all_estimates, na.rm = TRUE)
          mean_se <- colMeans(all_se, na.rm = TRUE)
          sd_estimates <- apply(all_estimates, 2, sd, na.rm = TRUE)
          se_sd_ratio <- mean_se / sd_estimates

          coverage <- numeric(n_coef)
          names(coverage) <- names(true_coef)
          for (coef_name in names(true_coef)) {
            in_ci <- (true_coef[coef_name] >= all_lower_ci[, coef_name]) &
                     (true_coef[coef_name] <= all_upper_ci[, coef_name])
            coverage[coef_name] <- mean(in_ci, na.rm = TRUE)
          }

          summary_table <- data.frame(
            Coefficient = names(true_coef),
            True_Value = true_coef,
            Mean_Estimate = mean_estimates,
            Mean_SE = mean_se,
            SD_Estimate = sd_estimates,
            SE_SD_Ratio = se_sd_ratio,
            Coverage = coverage,
            row.names = NULL,
            stringsAsFactors = FALSE
          )

          all_results[[scenario_name]] <- list(
            scenario = scenario_name,
            var_epsilon = var_epsilon,
            n = n,
            model_type = model_type,
            log_k = log_k,
            summary_table = summary_table,
            coverage = coverage,
            se_sd_ratio = se_sd_ratio,
            adjusted_count = adjusted_count,
            used_count = used_count
          )

          expect_true(
            all(coverage >= 0.85, na.rm = TRUE),
            info = sprintf(
              "Scenario %s: Coverage below 85%% for some coefficients.\n%s",
              scenario_name,
              paste(capture.output(print(summary_table[coverage < 0.85, ])), collapse = "\n")
            )
          )

          upper_se_sd <- 2.0
          if (n == 1e3 && var_epsilon == 0.9 && model_type == "weibull") {
            upper_se_sd <- 2.5
          }
          expect_true(
            all(se_sd_ratio >= 0.85 & se_sd_ratio <= upper_se_sd, na.rm = TRUE),
            info = sprintf(
              "Scenario %s: SE/SD ratio not within [0.85, %.1f] for some coefficients.\n%s",
              scenario_name,
              upper_se_sd,
              paste(capture.output(print(summary_table[se_sd_ratio < 0.85 | se_sd_ratio > upper_se_sd, ])), collapse = "\n")
            )
          )

          summary_file <- file.path(artifact_dir, paste0(scenario_name, "_summary.csv"))
          write.csv(summary_table, summary_file, row.names = FALSE)
          meta_file <- file.path(artifact_dir, paste0(scenario_name, "_meta.csv"))
          write.csv(
            data.frame(
              adjusted_count = adjusted_count,
              used_count = used_count,
              var_epsilon = var_epsilon,
              improvement_ratio = improvement_ratio,
              total_simulations = N_SIMULATIONS,
              stringsAsFactors = FALSE
            ),
            meta_file,
            row.names = FALSE
          )
        }
      }
    }
  }

  summary_all <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      n = x$n,
      var_epsilon = x$var_epsilon,
      model_type = x$model_type,
      log_k = x$log_k,
      Min_Coverage = min(x$coverage, na.rm = TRUE),
      Max_Coverage = max(x$coverage, na.rm = TRUE),
      Adjusted_Count = x$adjusted_count,
      Used_Count = x$used_count,
      stringsAsFactors = FALSE
    )
  }))

  summary_file <- file.path(artifact_dir, "summary.csv")
  write.csv(summary_all, summary_file, row.names = FALSE)
})
