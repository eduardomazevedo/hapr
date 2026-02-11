run_survival_tests <- function(test_type,
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
    true_coef <- make_true_coef(params$beta_g, params$beta_w)
    n_coef <- length(true_coef)

    theta <- normalize_theta(params$theta, scenario$var_v, scenario$var_epsilon)
    if (test_type == "point") {
      set.seed(scenario$seed)
      data <- if (scenario$model_type == "exponential") {
        mock_dataset_survival_exponential(
          n = scenario$n,
          var_v = scenario$var_v,
          var_epsilon = scenario$var_epsilon,
          beta_g = params$beta_g,
          beta_w = params$beta_w,
          theta = theta,
          censor_rate = params$censor_rate
        )
      } else {
        mock_dataset_survival_weibull(
          n = scenario$n,
          var_v = scenario$var_v,
          var_epsilon = scenario$var_epsilon,
          beta_g = params$beta_g,
          beta_w = params$beta_w,
          theta = theta,
          log_k = scenario$log_k,
          censor_rate = params$censor_rate
        )
      }

      start_beta <- rep(0, ncol(data$w) + 2)
      start_delta <- if (scenario$model_type == "weibull") c(log_k = 0) else numeric(0)

      mle_time <- system.time({
        mle_fit <- suppressWarnings(hapr_mle_survival(
          event_time = data$event_time,
          event_status = data$event_status,
          gc = data$gc,
          w = data$w,
          improvement_ratio = scenario$improvement_ratio,
          model_type = scenario$model_type,
          start_beta = start_beta,
          start_delta = start_delta,
          control = list(maxit = 150)
        ))
      })
      runtime_mle_ms <- mle_time[["elapsed"]] * 1000

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

      all_results[[scenario$name]] <- list(
        scenario = scenario$name,
        var_epsilon = scenario$var_epsilon,
        n = scenario$n,
        model_type = scenario$model_type,
        log_k = scenario$log_k,
        var_v_factor = scenario$var_v_factor,
        comparison_table = comparison_table,
        all_within_4se = all(within_4se, na.rm = TRUE)
      )

      expect_true(
        all(within_4se, na.rm = TRUE),
        info = sprintf(
          "Scenario %s: Some estimates are not within 4 SE of true values.\n%s",
          scenario$name,
          paste(capture.output(print(comparison_table[!within_4se & !is.na(within_4se), ])),
                collapse = "\n")
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

      adjusted_count <- 0
      used_count <- 0
      for (sim in 1:n_simulations) {
        set.seed(12345 +
                 scenario$idx_var_epsilon * 10000 +
                 scenario$idx_n * 1000 +
                 scenario$idx_model_type * 100 +
                 scenario$idx_var_v_factor * 10 +
                 round(scenario$log_k * 10) +
                 sim)

        data <- if (scenario$model_type == "exponential") {
          mock_dataset_survival_exponential(
            n = scenario$n,
            var_v = scenario$var_v,
            var_epsilon = scenario$var_epsilon,
            beta_g = params$beta_g,
            beta_w = params$beta_w,
            theta = theta,
            censor_rate = params$censor_rate
          )
        } else {
          mock_dataset_survival_weibull(
            n = scenario$n,
            var_v = scenario$var_v,
            var_epsilon = scenario$var_epsilon,
            beta_g = params$beta_g,
            beta_w = params$beta_w,
            theta = theta,
            log_k = scenario$log_k,
            censor_rate = params$censor_rate
          )
        }

        start_beta <- rep(0, ncol(data$w) + 2)
        start_delta <- if (scenario$model_type == "weibull") c(log_k = 0) else numeric(0)

        first_stage <- suppressWarnings(hapr_first_stage(
          y = data$event_time,
          gc = data$gc,
          w = data$w,
          model_type = "mle"
        ))
        max_improvement_ratio <- first_stage$stats$max_improvement_ratio
        if (scenario$n == 1e3 && scenario$improvement_ratio >= max_improvement_ratio) {
          adjusted_count <- adjusted_count + 1
          next
        }

        fit <- suppressWarnings(hapr_mle_survival(
          event_time = data$event_time,
          event_status = data$event_status,
          gc = data$gc,
          w = data$w,
          improvement_ratio = scenario$improvement_ratio,
          model_type = scenario$model_type,
          start_beta = start_beta,
          start_delta = start_delta,
          control = list(maxit = 150)
        ))

        ci_beta <- fit$ci_beta
        if (!is.null(ci_beta)) {
          aligned <- align_ci_beta(ci_beta, true_coef)
          all_estimates[sim, ] <- aligned$estimates
          all_se[sim, ] <- aligned$se
          all_lower_ci[sim, ] <- aligned$lower_ci
          all_upper_ci[sim, ] <- aligned$upper_ci
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

      all_results[[scenario$name]] <- list(
        scenario = scenario$name,
        var_epsilon = scenario$var_epsilon,
        n = scenario$n,
        model_type = scenario$model_type,
        log_k = scenario$log_k,
        var_v_factor = scenario$var_v_factor,
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
          scenario$name,
          paste(capture.output(print(summary_table[coverage < 0.85, ])), collapse = "\n")
        )
      )

      lower_se_sd <- 0.80
      upper_se_sd <- 3.0
      expect_true(
        all(se_sd_ratio >= lower_se_sd & se_sd_ratio <= upper_se_sd, na.rm = TRUE),
        info = sprintf(
          "Scenario %s: SE/SD ratio not within [%.2f, %.1f] for some coefficients.\n%s",
          scenario$name,
          lower_se_sd,
          upper_se_sd,
          paste(capture.output(print(summary_table[se_sd_ratio < lower_se_sd | se_sd_ratio > upper_se_sd, ])), collapse = "\n")
        )
      )

      write_summary_csv(artifact_dir, paste0(scenario$name, "_summary.csv"), summary_table)
      meta_file <- file.path(artifact_dir, paste0(scenario$name, "_meta.csv"))
      write.csv(
        data.frame(
          adjusted_count = adjusted_count,
          used_count = used_count,
          var_epsilon = scenario$var_epsilon,
          improvement_ratio = scenario$improvement_ratio,
          total_simulations = n_simulations,
          stringsAsFactors = FALSE
        ),
        meta_file,
        row.names = FALSE
      )
    }
  }

  all_results
}
