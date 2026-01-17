#' Test coverage intervals for hapr_second_stage across multiple scenarios
#' 
#' Runs 100 simulations per scenario and checks that coverage is above 85%
#' for all parameters. Tests the same scenarios as point estimate tests.
#' Produces artifact tables and histograms.
#' 
#' This test is slow and will only run if the environment variable
#' HAPR_RUN_COVERAGE_TESTS is set to "true" or "TRUE".

test_that("Coverage intervals are above 85% for all scenarios", {
  # Check if coverage tests should run
  run_coverage <- Sys.getenv("HAPR_RUN_COVERAGE_TESTS", unset = "false")
  if (!tolower(run_coverage) %in% c("true", "1", "yes")) {
    skip("Coverage tests skipped. Set HAPR_RUN_COVERAGE_TESTS=true to run.")
  }
  
  # Source mock dataset functions (using path relative to project root)
  # testthat runs from project root
  mock_dataset_path <- file.path("dev", "mock_dataset.R")
  if (file.exists(mock_dataset_path)) {
    source(mock_dataset_path)
  } else {
    # Try alternative path if running from testthat directory
    mock_dataset_path <- file.path("..", "..", "dev", "mock_dataset.R")
    if (file.exists(mock_dataset_path)) {
      source(mock_dataset_path)
    } else {
      stop("Cannot find dev/mock_dataset.R")
    }
  }
  
  # True coefficients (same as in stage_2_lm.R)
  BETA_G <- 1.42  # Effect of future PRS (gf) on outcome
  BETA_W <- c(0.1, 0.17, 0.27, -0.27)  # Intercept + 3 covariate effects
  THETA <- c(0.0, 0.1, -0.2, 0.3)  # Intercept + 3 covariate effects for gc ~ w
  VAR_Y <- 1.0  # Error variance for outcome (only for lm)
  
  # Test scenarios
  VAR_EPSILON_VALUES <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  N_VALUES <- c(1e3, 1e4, 1e5)
  MODEL_TYPES <- c("lm", "probit")
  N_SIMULATIONS <- 100
  
  # Create directory for artifacts if it doesn't exist
  artifact_dir <- testthat::test_path("_artifacts", "coverage")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }
  
  # Store results for all scenarios
  all_results <- list()
  
  for (var_epsilon in VAR_EPSILON_VALUES) {
    for (n in N_VALUES) {
      for (model_type in MODEL_TYPES) {
        scenario_name <- sprintf("%s_n%d_ve%.1f", model_type, as.integer(n), var_epsilon)
        
        # Skip probit with n=1000 and var_epsilon=0.9 - sample size too small
        if (model_type == "probit" && n == 1e3 && var_epsilon == 0.9) {
          next
        }
        
        # Calculate var_v ensuring that 1 - var_v - var_epsilon >= 0
        var_v <- (1 - var_epsilon) * 0.5
        improvement_ratio <- 1 / (1 - var_epsilon)
        
        # Construct true coefficients vector
        true_coef <- c(BETA_G, BETA_W)
        names(true_coef) <- c("gf", "(Intercept)", "w1", "w2", "w3")
        n_coef <- length(true_coef)
        
        # Storage for simulation results
        all_estimates <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
        colnames(all_estimates) <- names(true_coef)
        all_se <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
        colnames(all_se) <- names(true_coef)
        all_lower_ci <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
        colnames(all_lower_ci) <- names(true_coef)
        all_upper_ci <- matrix(NA, nrow = N_SIMULATIONS, ncol = n_coef)
        colnames(all_upper_ci) <- names(true_coef)
        
        # Run simulations
        for (sim in 1:N_SIMULATIONS) {
          # Set seed for reproducibility (different seed per simulation)
          set.seed(12345 + 
                   which(VAR_EPSILON_VALUES == var_epsilon) * 10000 + 
                   which(N_VALUES == n) * 1000 + 
                   which(MODEL_TYPES == model_type) * 100 + 
                   sim)
          
          # Create dataset
          if (model_type == "lm") {
            data <- mock_dataset_lm(
              n = n,
              var_v = var_v,
              var_epsilon = var_epsilon,
              beta_g = BETA_G,
              beta_w = BETA_W,
              theta = THETA,
              var_y = VAR_Y
            )
          } else {  # probit
            data <- mock_dataset_probit(
              n = n,
              var_v = var_v,
              var_epsilon = var_epsilon,
              beta_g = BETA_G,
              beta_w = BETA_W,
              theta = THETA
            )
          }
          
          w <- data$w
          gc <- data$gc
          y <- data$y
          
          # Run hapr stages
          first_stage_fit <- hapr_first_stage(
            y = y,
            gc = gc,
            w = w,
            model_type = model_type
          )
          
          # Try second stage with error handling
          tryCatch({
            second_stage_fit <- hapr_second_stage(
              first_stage = first_stage_fit,
              improvement_ratio = improvement_ratio
            )
            
            # Extract estimates and standard errors
            ci_beta <- second_stage_fit$ci_beta
            coef_names <- rownames(ci_beta)
            estimates <- setNames(ci_beta$Estimate, coef_names)
            se <- setNames(ci_beta$Std.Error, coef_names)
            lower_ci <- setNames(ci_beta$Lower, coef_names)
            upper_ci <- setNames(ci_beta$Upper, coef_names)
            
            # Store results (ordered by true_coef names)
            all_estimates[sim, ] <- estimates[names(true_coef)]
            all_se[sim, ] <- se[names(true_coef)]
            all_lower_ci[sim, ] <- lower_ci[names(true_coef)]
            all_upper_ci[sim, ] <- upper_ci[names(true_coef)]
          }, error = function(e) {
            # Report the exact scenario and simulation that failed
            error_msg <- sprintf(
              "Error in scenario %s, simulation %d:\n  model_type=%s, n=%d, var_epsilon=%.2f, var_v=%.4f, improvement_ratio=%.4f\n  Error: %s",
              scenario_name, sim, model_type, n, var_epsilon, var_v, improvement_ratio, e$message
            )
            # Also print first stage gamma_gc for debugging
            gamma_gc <- first_stage_fit$parameters$gamma[1]
            cat("\n", error_msg, "\n")
            cat(sprintf("  First stage gamma_gc: %.6f\n", gamma_gc))
            cat(sprintf("  First stage var_v_plus_var_epsilon: %.6f\n", first_stage_fit$parameters$var_v_plus_var_epsilon))
            stop(error_msg)
          })
        }
        
        # Calculate statistics for each coefficient
        mean_estimates <- colMeans(all_estimates, na.rm = TRUE)
        mean_se <- colMeans(all_se, na.rm = TRUE)
        sd_estimates <- apply(all_estimates, 2, sd, na.rm = TRUE)
        se_sd_ratio <- mean_se / sd_estimates
        
        # Calculate coverage (proportion of simulations where true value is in CI)
        coverage <- numeric(n_coef)
        names(coverage) <- names(true_coef)
        for (coef_name in names(true_coef)) {
          in_ci <- (true_coef[coef_name] >= all_lower_ci[, coef_name]) & 
                   (true_coef[coef_name] <= all_upper_ci[, coef_name])
          coverage[coef_name] <- mean(in_ci, na.rm = TRUE)
        }
        
        # Create summary table
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
        
        # Calculate mean CI bounds for histograms
        mean_lower_ci <- colMeans(all_lower_ci, na.rm = TRUE)
        mean_upper_ci <- colMeans(all_upper_ci, na.rm = TRUE)
        
        # Store results
        all_results[[scenario_name]] <- list(
          scenario = scenario_name,
          var_epsilon = var_epsilon,
          n = n,
          model_type = model_type,
          summary_table = summary_table,
          all_estimates = all_estimates,
          mean_lower_ci = mean_lower_ci,
          mean_upper_ci = mean_upper_ci,
          coverage = coverage,
          se_sd_ratio = se_sd_ratio
        )
        
        # Test assertions
        # 1. Coverage should be above 85% for all coefficients
        expect_true(
          all(coverage >= 0.85, na.rm = TRUE),
          info = sprintf(
            "Scenario %s: Coverage below 85%% for some coefficients.\n%s",
            scenario_name,
            paste(capture.output(print(summary_table[coverage < 0.85, ])), collapse = "\n")
          )
        )
        
        # 2. SE/SD ratio should be between 0.85 and 2.0
        expect_true(
          all(se_sd_ratio >= 0.85 & se_sd_ratio <= 2.0, na.rm = TRUE),
          info = sprintf(
            "Scenario %s: SE/SD ratio not within [0.85, 2.0] for some coefficients.\n%s",
            scenario_name,
            paste(capture.output(print(summary_table[se_sd_ratio < 0.85 | se_sd_ratio > 2.0, ])), collapse = "\n")
          )
        )
        
        # Save summary table
        summary_file <- file.path(artifact_dir, paste0(scenario_name, "_summary.csv"))
        write.csv(summary_table, summary_file, row.names = FALSE)
        
        # Generate histograms for each coefficient
        for (coef_name in names(true_coef)) {
          estimates_vec <- all_estimates[, coef_name]
          estimates_vec <- estimates_vec[!is.na(estimates_vec)]
          
          if (length(estimates_vec) > 0) {
            # Create histogram
            png_file <- file.path(artifact_dir, paste0(scenario_name, "_", coef_name, "_hist.png"))
            png(png_file, width = 800, height = 600, res = 100)
            
            # Create histogram and store result to get max count
            h <- hist(estimates_vec, 
                     main = sprintf("%s: %s\nTrue=%.3f, Mean=%.3f, Coverage=%.1f%%",
                                   scenario_name, coef_name,
                                   true_coef[coef_name], mean_estimates[coef_name],
                                   coverage[coef_name] * 100),
                     xlab = "Estimate",
                     ylab = "Frequency",
                     breaks = 20,
                     col = "lightblue",
                     border = "black")
            
            # Add vertical line for true value
            abline(v = true_coef[coef_name], col = "red", lwd = 2, lty = 2)
            
            # Add vertical line for mean estimate
            abline(v = mean_estimates[coef_name], col = "blue", lwd = 2, lty = 2)
            
            # Add bars for mean 95% CI (horizontal arrow at 90% of max count)
            max_count <- max(h$counts)
            arrows(mean_lower_ci[coef_name], max_count * 0.9,
                   mean_upper_ci[coef_name], max_count * 0.9,
                   code = 3, angle = 90, length = 0.1, col = "green", lwd = 3)
            
            # Add legend
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
  }
  
  # Create overall summary table
  overall_summary <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      Model_Type = x$model_type,
      n = x$n,
      var_epsilon = x$var_epsilon,
      Min_Coverage = min(x$coverage, na.rm = TRUE),
      Max_Coverage = max(x$coverage, na.rm = TRUE),
      Mean_Coverage = mean(x$coverage, na.rm = TRUE),
      Min_SE_SD_Ratio = min(x$se_sd_ratio, na.rm = TRUE),
      Max_SE_SD_Ratio = max(x$se_sd_ratio, na.rm = TRUE),
      Mean_SE_SD_Ratio = mean(x$se_sd_ratio, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  
  # Save overall summary
  overall_summary_file <- file.path(artifact_dir, "overall_summary.csv")
  write.csv(overall_summary, overall_summary_file, row.names = FALSE)
  
  # Final test: all scenarios should have coverage >= 85% and SE/SD ratio in [0.85, 2.0]
  all_coverage_ok <- all(sapply(all_results, function(x) all(x$coverage >= 0.85, na.rm = TRUE)))
  all_ratio_ok <- all(sapply(all_results, function(x) all(x$se_sd_ratio >= 0.85 & x$se_sd_ratio <= 2.0, na.rm = TRUE)))
  
  expect_true(
    all_coverage_ok && all_ratio_ok,
    info = sprintf(
      "Not all scenarios passed. Overall summary:\n%s",
      paste(capture.output(print(overall_summary)), collapse = "\n")
    )
  )
})
