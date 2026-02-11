test_that("Aggregate coverage artifacts", {
  run_slow <- Sys.getenv("RUN_SLOW_TESTS", unset = "false")
  if (!tolower(run_slow) %in% c("true", "1", "yes")) {
    skip("Aggregation skipped. Set RUN_SLOW_TESTS=true to run.")
  }

  pkg_root <- testthat::test_path("..", "..")
  dev_script <- file.path(pkg_root, "dev", "aggregate_artifacts.R")
  if (!file.exists(dev_script)) {
    skip("Aggregation script not found.")
  }

  source(dev_script)
  artifact_root <- testthat::test_path("_artifacts")
  output_file <- file.path(artifact_root, "coverage_summary.csv")
  aggregate_coverage_artifacts(artifact_root = artifact_root, output_file = output_file)

  expect_true(file.exists(output_file))
})
