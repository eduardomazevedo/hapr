#' Test point estimates for hapr_second_stage across multiple scenarios
#' 
#' Tests that point estimates are within 3 standard errors of true coefficients
#' for various combinations of var_epsilon, n, and model_type.
#' Produces artifact tables comparing true coefficients, estimates, SEs, and CIs.

test_that("Point estimates are within 3 SE of true coefficients for all scenarios", {
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
  
  # Create directory for artifacts if it doesn't exist
  # Note: testthat runs from project root, so paths are relative to that
  artifact_dir <- testthat::test_path("_artifacts", "point_estimates")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }
  # Store results for all scenarios
  all_results <- list()
  
  for (var_epsilon in VAR_EPSILON_VALUES) {
    for (n in N_VALUES) {
      for (model_type in MODEL_TYPES) {
        scenario_name <- sprintf("%s_n%d_ve%.1f", model_type, as.integer(n), var_epsilon)
        
        # Set seed for reproducibility (different seed per scenario)
        set.seed(123 + which(VAR_EPSILON_VALUES == var_epsilon) * 100 + 
                        which(N_VALUES == n) * 10 + 
                        which(MODEL_TYPES == model_type))
        
        # Calculate var_v ensuring that 1 - var_v - var_epsilon >= 0
        # This is required for the normalization in mock_dataset_lm
        # Use a formula that keeps var_v + var_epsilon < 1
        # For example: var_v = (1 - var_epsilon) * 0.5 ensures var_v + var_epsilon = 0.5 + 0.5*var_epsilon < 1
        var_v <- (1 - var_epsilon) * 0.5
        improvement_ratio <- 1 / (1 - var_epsilon)
        
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
        two_stage_time <- system.time({
          first_stage_fit <- hapr_first_stage(
            y = y,
            gc = gc,
            w = w,
            model_type = model_type
          )
          
          second_stage_fit <- hapr_second_stage(
            first_stage = first_stage_fit,
            improvement_ratio = improvement_ratio
          )
        })
        runtime_two_stage_ms <- two_stage_time[["elapsed"]] * 1000
        
        # Extract estimates and standard errors
        ci_beta <- second_stage_fit$ci_beta
        # ci_beta is a data frame with row names as coefficient names
        coef_names <- rownames(ci_beta)
        estimates <- setNames(ci_beta$Estimate, coef_names)
        se <- setNames(ci_beta$Std.Error, coef_names)
        lower_ci <- setNames(ci_beta$Lower, coef_names)
        upper_ci <- setNames(ci_beta$Upper, coef_names)
        
        # Construct true coefficients vector (matching the order in ci_beta)
        # Beta coefficients are named: "gf", "(Intercept)", "w1", "w2", "w3"
        true_coef <- c(BETA_G, BETA_W)
        names(true_coef) <- c("gf", "(Intercept)", "w1", "w2", "w3")
        
        # Ensure estimates are in the same order as true_coef
        estimates_ordered <- estimates[names(true_coef)]
        se_ordered <- se[names(true_coef)]
        lower_ci_ordered <- lower_ci[names(true_coef)]
        upper_ci_ordered <- upper_ci[names(true_coef)]
        
        # Check that all estimates are within 3 SE of true values
        differences <- abs(estimates_ordered - true_coef)
        within_3se <- differences <= 3 * se_ordered
        
        # Create comparison table
        comparison_table <- data.frame(
          Coefficient = names(true_coef),
          True_Value = true_coef,
          Point_Estimate = estimates_ordered,
          Std_Error = se_ordered,
          Lower_CI = lower_ci_ordered,
          Upper_CI = upper_ci_ordered,
          Difference = differences,
          Within_3SE = within_3se,
          Runtime_Two_Stage_ms = rep(runtime_two_stage_ms, length(true_coef)),
          row.names = NULL,
          stringsAsFactors = FALSE
        )
        
        # Store results
        all_results[[scenario_name]] <- list(
          scenario = scenario_name,
          var_epsilon = var_epsilon,
          n = n,
          model_type = model_type,
          comparison_table = comparison_table,
          all_within_3se = all(within_3se, na.rm = TRUE)
        )
        
        # Test assertion: all estimates should be within 3 SE
        expect_true(
          all(within_3se, na.rm = TRUE),
          info = sprintf(
            "Scenario %s: Some estimates are not within 3 SE of true values.\n%s",
            scenario_name,
            paste(capture.output(print(comparison_table[!within_3se & !is.na(within_3se), ])), collapse = "\n")
          )
        )
        
        # Save artifact table
        artifact_file <- file.path(artifact_dir, paste0(scenario_name, ".csv"))
        write.csv(comparison_table, artifact_file, row.names = FALSE)
      }
    }
  }
  
  # Create summary table with all scenarios
  summary_table <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      Model_Type = x$model_type,
      n = x$n,
      var_epsilon = x$var_epsilon,
      All_Within_3SE = x$all_within_3se,
      stringsAsFactors = FALSE
    )
  }))
  
  # Save summary artifact
  summary_file <- file.path(artifact_dir, "summary.csv")
  write.csv(summary_table, summary_file, row.names = FALSE)
  
  # Test that all scenarios passed
  expect_true(
    all(sapply(all_results, function(x) x$all_within_3se), na.rm = TRUE),
    info = sprintf(
      "Not all scenarios passed. Summary:\n%s",
      paste(capture.output(print(summary_table)), collapse = "\n")
    )
  )
})
