ensure_artifact_dir <- function(...) {
  artifact_dir <- testthat::test_path("_artifacts", ...)
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }
  artifact_dir
}

write_scenario_csv <- function(artifact_dir, scenario_name, data) {
  artifact_file <- file.path(artifact_dir, paste0(scenario_name, ".csv"))
  write.csv(data, artifact_file, row.names = FALSE)
}

write_summary_csv <- function(artifact_dir, filename, data) {
  summary_file <- file.path(artifact_dir, filename)
  write.csv(data, summary_file, row.names = FALSE)
}
