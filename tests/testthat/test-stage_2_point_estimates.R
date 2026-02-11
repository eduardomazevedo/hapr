#' Test point estimates for hapr_second_stage across multiple scenarios
#' 
#' Tests that point estimates are within 3 standard errors of true coefficients
#' for various combinations of var_epsilon, n, and model_type.
#' Produces artifact tables comparing true coefficients, estimates, SEs, and CIs.

test_that("Point estimates are within 3 SE of true coefficients for all scenarios", {
  params <- stage2_params_default()
  scenarios <- stage2_scenarios(is_slow_enabled(), include_large_n = TRUE)
  artifact_dir <- ensure_artifact_dir("point_estimates")

  all_results <- run_stage2_tests(
    test_type = "point",
    params = params,
    scenarios = scenarios,
    artifact_dir = artifact_dir
  )

  summary_table <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      Model_Type = x$model_type,
      n = x$n,
      var_epsilon = x$var_epsilon,
      All_Within_3SE = x$all_within_3se,
      stringsAsFactors = FALSE
    )
  }))

  write_summary_csv(artifact_dir, "summary.csv", summary_table)

  expect_true(
    all(sapply(all_results, function(x) x$all_within_3se), na.rm = TRUE),
    info = sprintf(
      "Not all scenarios passed. Summary:\n%s",
      paste(capture.output(print(summary_table)), collapse = "\n")
    )
  )
})
