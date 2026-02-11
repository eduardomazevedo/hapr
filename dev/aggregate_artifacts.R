aggregate_coverage_artifacts <- function(
  artifact_root = file.path(getwd(), "tests", "testthat", "_artifacts"),
  output_file = file.path(artifact_root, "coverage_summary.csv")
) {
  coverage_dirs <- c("coverage", "coverage_survival")
  csv_paths <- unlist(lapply(coverage_dirs, function(dir_name) {
    dir_path <- file.path(artifact_root, dir_name)
    if (!dir.exists(dir_path)) {
      return(character(0))
    }
    list.files(dir_path, pattern = "_summary\\.csv$", full.names = TRUE)
  }), use.names = FALSE)

  if (!length(csv_paths)) {
    message("No summary CSVs found under ", artifact_root)
    return(invisible(NULL))
  }

  csv_paths <- csv_paths[!grepl("(/|\\\\)(summary|overall_summary)\\.csv$", csv_paths)]
  if (!length(csv_paths)) {
    message("No scenario summary CSVs found under ", artifact_root)
    return(invisible(NULL))
  }

  parse_scenario <- function(scenario) {
    result <- list(
      model_type = NA_character_,
      n = NA_integer_,
      var_epsilon = NA_real_,
      var_v_factor = NA_real_,
      log_k = NA_real_
    )

    coverage_match <- regexec("^(lm|probit)_n([0-9]+)_ve([0-9.]+)_vv([0-9.]+)$", scenario)
    coverage_parts <- regmatches(scenario, coverage_match)[[1]]
    if (length(coverage_parts) > 0) {
      result$model_type <- coverage_parts[2]
      result$n <- as.integer(coverage_parts[3])
      result$var_epsilon <- as.numeric(coverage_parts[4])
      result$var_v_factor <- as.numeric(coverage_parts[5])
      return(result)
    }

    survival_match <- regexec("^(exp|wei)_n([0-9]+)_ve([0-9.]+)_vv([0-9.]+)_lk(-?[0-9.]+)$", scenario)
    survival_parts <- regmatches(scenario, survival_match)[[1]]
    if (length(survival_parts) > 0) {
      result$model_type <- survival_parts[2]
      result$n <- as.integer(survival_parts[3])
      result$var_epsilon <- as.numeric(survival_parts[4])
      result$var_v_factor <- as.numeric(survival_parts[5])
      result$log_k <- as.numeric(survival_parts[6])
    }

    result
  }

  combined <- do.call(rbind, lapply(csv_paths, function(path) {
    scenario <- sub("_summary\\.csv$", "", basename(path))
    suite <- basename(dirname(path))
    summary_table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    meta <- parse_scenario(scenario)

    meta_table <- data.frame(
      suite = suite,
      scenario = scenario,
      model_type = meta$model_type,
      n = meta$n,
      var_epsilon = meta$var_epsilon,
      var_v_factor = meta$var_v_factor,
      log_k = meta$log_k,
      stringsAsFactors = FALSE
    )

    # Ensure numeric columns stay numeric
    for (col in names(summary_table)) {
      if (grepl("^(True_Value|Mean_Estimate|Mean_SE|SD_Estimate|SE_SD_Ratio|Coverage)$", col)) {
        summary_table[[col]] <- as.numeric(summary_table[[col]])
      }
    }

    meta_table <- meta_table[rep(1, nrow(summary_table)), , drop = FALSE]
    rownames(meta_table) <- NULL
    cbind(meta_table, summary_table, row.names = NULL)
  }))

  rownames(combined) <- NULL
  write.csv(combined, output_file, row.names = FALSE)
  invisible(combined)
}

if (sys.nframe() == 0) {
  aggregate_coverage_artifacts()
}
