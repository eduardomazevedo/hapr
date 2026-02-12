#' Coverage tests for survival MLE
#'
#' Runs 100 (or 1e3 if RUN_GIANT_TESTS is set) simulations per scenario and
#' checks that coverage is above 85% for all beta parameters. This test is
#' slow and only runs when RUN_SLOW_TESTS is set. Set RUN_GIANT_TESTS=true
#' for n = 1e5 and 1e3 simulations.

test_that("Survival MLE coverage intervals are above 85%", {
  if (!is_slow_enabled()) {
    skip("Coverage tests skipped. Set RUN_SLOW_TESTS=true or RUN_GIANT_TESTS=true to run.")
  }

  params <- survival_params_default()
  scenarios <- survival_scenarios(run_slow = TRUE, log_k_values = c(0), include_large_n = is_giant_enabled())
  artifact_dir <- ensure_artifact_dir("coverage_survival")

  all_results <- run_survival_tests(
    test_type = "coverage",
    params = params,
    scenarios = scenarios,
    artifact_dir = artifact_dir,
    n_simulations = if (is_giant_enabled()) 1e3 else 100
  )

  summary_all <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      n = x$n,
      var_epsilon = x$var_epsilon,
      model_type = x$model_type,
      log_k = x$log_k,
      var_v_factor = x$var_v_factor,
      Min_Coverage = min(x$coverage, na.rm = TRUE),
      Max_Coverage = max(x$coverage, na.rm = TRUE),
      Adjusted_Count = x$adjusted_count,
      Used_Count = x$used_count,
      stringsAsFactors = FALSE
    )
  }))

  write_summary_csv(artifact_dir, "summary.csv", summary_all)
})
