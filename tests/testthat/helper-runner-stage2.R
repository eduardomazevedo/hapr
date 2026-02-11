run_stage2_tests <- function(test_type,
                             params,
                             scenarios,
                             artifact_dir,
                             n_simulations = NULL) {
  if (!test_type %in% c("point", "coverage")) {
    stop("test_type must be one of: 'point', 'coverage'")
  }

  if (test_type == "coverage" && is.null(n_simulations)) {
    stop("n_simulations must be set for coverage tests")
  }

  all_results <- list()

  for (scenario in scenarios) {
    # Skip probit with n=1000 and var_epsilon=0.9 - sample size too small
    if (scenario$model_type == "probit" && scenario$n == 1e3 && scenario$var_epsilon == 0.9) {
      next
    }

    true_coef <- make_true_coef(params$beta_g, params$beta_w)
    n_coef <- length(true_coef)

    theta <- normalize_theta(params$theta, scenario$var_v, scenario$var_epsilon)
    if (test_type == "point") {
      set.seed(scenario$seed)
      data <- if (scenario$model_type == "lm") {
        mock_dataset_lm(
          n = scenario$n,
          var_v = scenario$var_v,
          var_epsilon = scenario$var_epsilon,
          beta_g = params$beta_g,
          beta_w = params$beta_w,
          theta = theta,
          var_y = params$var_y
        )
      } else {
        mock_dataset_probit(
          n = scenario$n,
          var_v = scenario$var_v,
          var_epsilon = scenario$var_epsilon,
          beta_g = params$beta_g,
          beta_w = params$beta_w,
          theta = theta
        )
      }

      two_stage_time <- system.time({
        first_stage_fit <- suppressWarnings(hapr_first_stage(
          y = data$y,
          gc = data$gc,
          w = data$w,
          model_type = scenario$model_type
        ))
        second_stage_fit <- suppressWarnings(hapr_second_stage(
          first_stage = first_stage_fit,
          improvement_ratio = scenario$improvement_ratio
        ))
      })
      runtime_two_stage_ms <- two_stage_time[["elapsed"]] * 1000

      aligned <- align_ci_beta(second_stage_fit$ci_beta, true_coef)
      differences <- abs(aligned$estimates - true_coef)
      within_3se <- differences <= 3 * aligned$se

      comparison_table <- data.frame(
        Coefficient = names(true_coef),
        True_Value = true_coef,
        Point_Estimate = aligned$estimates,
        Std_Error = aligned$se,
        Lower_CI = aligned$lower_ci,
        Upper_CI = aligned$upper_ci,
        Difference = differences,
        Within_3SE = within_3se,
        Runtime_Two_Stage_ms = rep(runtime_two_stage_ms, length(true_coef)),
        row.names = NULL,
        stringsAsFactors = FALSE
      )

      all_results[[scenario$name]] <- list(
        scenario = scenario$name,
        var_epsilon = scenario$var_epsilon,
        n = scenario$n,
        model_type = scenario$model_type,
        var_v_factor = scenario$var_v_factor,
        comparison_table = comparison_table,
        all_within_3se = all(within_3se, na.rm = TRUE)
      )

      expect_true(
        all(within_3se, na.rm = TRUE),
        info = sprintf(
          "Scenario %s: Some estimates are not within 3 SE of true values.\n%s",
          scenario$name,
          paste(capture.output(print(comparison_table[!within_3se & !is.na(within_3se), ])), collapse = "\n")
        )
      )

      write_scenario_csv(artifact_dir, scenario$name, comparison_table)
    } else {
      all_estimates <- matrix(NA, nrow = n_simulations, ncol = n_coef)
      colnames(all_estimates) <- names(true_coef)
      all_se <- matrix(NA, nrow = n_simulations, ncol = n_coef)
      colnames(all_se) <- names(true_coef)
      all_lower_ci <- matrix(NA, nrow = n_simulations, ncol = n_coef)
      colnames(all_lower_ci) <- names(true_coef)
      all_upper_ci <- matrix(NA, nrow = n_simulations, ncol = n_coef)
      colnames(all_upper_ci) <- names(true_coef)

      for (sim in 1:n_simulations) {
        set.seed(12345 +
                   scenario$idx_var_epsilon * 10000 +
                   scenario$idx_n * 1000 +
                   scenario$idx_model_type * 100 +
                   scenario$idx_var_v_factor * 10 +
                   sim)

        data <- if (scenario$model_type == "lm") {
          mock_dataset_lm(
            n = scenario$n,
            var_v = scenario$var_v,
            var_epsilon = scenario$var_epsilon,
            beta_g = params$beta_g,
            beta_w = params$beta_w,
            theta = theta,
            var_y = params$var_y
          )
        } else {
          mock_dataset_probit(
            n = scenario$n,
            var_v = scenario$var_v,
            var_epsilon = scenario$var_epsilon,
            beta_g = params$beta_g,
            beta_w = params$beta_w,
            theta = theta
          )
        }

        first_stage_fit <- suppressWarnings(hapr_first_stage(
          y = data$y,
          gc = data$gc,
          w = data$w,
          model_type = scenario$model_type
        ))

        tryCatch({
          second_stage_fit <- suppressWarnings(hapr_second_stage(
            first_stage = first_stage_fit,
            improvement_ratio = scenario$improvement_ratio
          ))

          aligned <- align_ci_beta(second_stage_fit$ci_beta, true_coef)
          all_estimates[sim, ] <- aligned$estimates
          all_se[sim, ] <- aligned$se
          all_lower_ci[sim, ] <- aligned$lower_ci
          all_upper_ci[sim, ] <- aligned$upper_ci
        }, error = function(e) {
          error_msg <- sprintf(
            "Error in scenario %s, simulation %d:\n  model_type=%s, n=%d, var_epsilon=%.2f, var_v=%.4f, improvement_ratio=%.4f\n  Error: %s",
            scenario$name,
            sim,
            scenario$model_type,
            scenario$n,
            scenario$var_epsilon,
            scenario$var_v,
            scenario$improvement_ratio,
            e$message
          )
          gamma_gc <- first_stage_fit$parameters$gamma[1]
          cat("\n", error_msg, "\n")
          cat(sprintf("  First stage gamma_gc: %.6f\n", gamma_gc))
          cat(sprintf("  First stage var_v_plus_var_epsilon: %.6f\n", first_stage_fit$parameters$var_v_plus_var_epsilon))
          stop(error_msg)
        })
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

      mean_lower_ci <- colMeans(all_lower_ci, na.rm = TRUE)
      mean_upper_ci <- colMeans(all_upper_ci, na.rm = TRUE)

      all_results[[scenario$name]] <- list(
        scenario = scenario$name,
        var_epsilon = scenario$var_epsilon,
        n = scenario$n,
        model_type = scenario$model_type,
        var_v_factor = scenario$var_v_factor,
        summary_table = summary_table,
        all_estimates = all_estimates,
        mean_lower_ci = mean_lower_ci,
        mean_upper_ci = mean_upper_ci,
        coverage = coverage,
        se_sd_ratio = se_sd_ratio
      )

      expect_true(
        all(coverage >= 0.85, na.rm = TRUE),
        info = sprintf(
          "Scenario %s: Coverage below 85%% for some coefficients.\n%s",
          scenario$name,
          paste(capture.output(print(summary_table[coverage < 0.85, ])), collapse = "\n")
        )
      )

      lower_se_sd <- 0.80
      upper_se_sd <- 2.30
      expect_true(
        all(se_sd_ratio >= lower_se_sd & se_sd_ratio <= upper_se_sd, na.rm = TRUE),
        info = sprintf(
          "Scenario %s: SE/SD ratio not within [%.2f, %.2f] for some coefficients.\n%s",
          scenario$name,
          lower_se_sd,
          upper_se_sd,
          paste(capture.output(print(summary_table[se_sd_ratio < lower_se_sd | se_sd_ratio > upper_se_sd, ])), collapse = "\n")
        )
      )

      write_summary_csv(artifact_dir, paste0(scenario$name, "_summary.csv"), summary_table)

      for (coef_name in names(true_coef)) {
        estimates_vec <- all_estimates[, coef_name]
        estimates_vec <- estimates_vec[!is.na(estimates_vec)]

        if (length(estimates_vec) > 0) {
          png_file <- file.path(artifact_dir, paste0(scenario$name, "_", coef_name, "_hist.png"))
          png(png_file, width = 800, height = 600, res = 100)

          h <- hist(estimates_vec,
                    main = sprintf("%s: %s\nTrue=%.3f, Mean=%.3f, Coverage=%.1f%%",
                                   scenario$name, coef_name,
                                   true_coef[coef_name], mean_estimates[coef_name],
                                   coverage[coef_name] * 100),
                    xlab = "Estimate",
                    ylab = "Frequency",
                    breaks = 20,
                    col = "lightblue",
                    border = "black")

          abline(v = true_coef[coef_name], col = "red", lwd = 2, lty = 2)
          abline(v = mean_estimates[coef_name], col = "blue", lwd = 2, lty = 2)

          max_count <- max(h$counts)
          arrows(mean_lower_ci[coef_name], max_count * 0.9,
                 mean_upper_ci[coef_name], max_count * 0.9,
                 code = 3, angle = 90, length = 0.1, col = "green", lwd = 3)

          legend("topright",
                 legend = c(sprintf("True: %.3f", true_coef[coef_name]),
                            sprintf("Mean: %.3f", mean_estimates[coef_name]),
                            sprintf("Mean 95%% CI: [%.3f, %.3f]",
                                    mean_lower_ci[coef_name], mean_upper_ci[coef_name])),
                 col = c("red", "blue", "green"),
                 lty = c(2, 2, 1),
                 lwd = c(2, 2, 3))

          dev.off()
        }
      }
    }
  }

  all_results
}
