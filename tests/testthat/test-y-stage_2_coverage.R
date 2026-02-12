#' Test coverage intervals for hapr_second_stage across multiple scenarios
#' 
#' Runs 1e3 simulations per scenario and checks that coverage is above 85%
#' for all parameters. Tests the same scenarios as point estimate tests.
#' Produces artifact tables.
#' 
#' This test is slow and will only run if the environment variable
#' RUN_SLOW_TESTS is set to "true" or "TRUE".

test_that("Coverage intervals are above 85% for all scenarios", {
  # Check if coverage tests should run
  if (!is_slow_enabled()) {
    skip("Coverage tests skipped. Set RUN_SLOW_TESTS=true to run.")
  }

  params <- stage2_params_default()
  scenarios <- stage2_scenarios(run_slow = TRUE, include_large_n = TRUE)
  artifact_dir <- ensure_artifact_dir("coverage")

  all_results <- run_stage2_tests(
    test_type = "coverage",
    params = params,
    scenarios = scenarios,
    artifact_dir = artifact_dir,
    n_simulations = 1e3
  )

  overall_summary <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      Model_Type = x$model_type,
      n = x$n,
      var_epsilon = x$var_epsilon,
      var_v_factor = x$var_v_factor,
      Min_Coverage = min(x$coverage, na.rm = TRUE),
      Max_Coverage = max(x$coverage, na.rm = TRUE),
      Mean_Coverage = mean(x$coverage, na.rm = TRUE),
      Min_SE_SD_Ratio = min(x$se_sd_ratio, na.rm = TRUE),
      Max_SE_SD_Ratio = max(x$se_sd_ratio, na.rm = TRUE),
      Mean_SE_SD_Ratio = mean(x$se_sd_ratio, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  write_summary_csv(artifact_dir, "overall_summary.csv", overall_summary)

  lower_se_sd <- 0.80
  upper_se_sd <- 2.30
  all_coverage_ok <- all(sapply(all_results, function(x) all(x$coverage >= 0.85, na.rm = TRUE)))
  all_ratio_ok <- all(sapply(
    all_results,
    function(x) all(x$se_sd_ratio >= lower_se_sd & x$se_sd_ratio <= upper_se_sd, na.rm = TRUE)
  ))

  expect_true(
    all_coverage_ok && all_ratio_ok,
    info = sprintf(
      "Not all scenarios passed. Overall summary:\n%s",
      paste(capture.output(print(overall_summary)), collapse = "\n")
    )
  )
})
