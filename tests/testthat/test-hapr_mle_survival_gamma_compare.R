#' Compare gamma-based survival MLE against beta-based MLE
#'
#' Ensures point estimates match closely, records runtimes, and checks
#' that beta SEs from gamma-based delta method are larger.

test_that("Gamma-based survival MLE matches beta-based MLE", {
  set.seed(321)

  n <- 100000
  var_epsilon <- 0.7
  var_v <- (1 - var_epsilon) * 0.4
  improvement_ratio <- 1 / (1 - var_epsilon)

  beta_g <- 0.6
  beta_w <- c(0.1, -0.2, 0.15, 0.05)
  theta <- c(0.0, 0.1, -0.25, 0.2)
  censor_rate <- 0.2

  model_types <- c("exponential", "weibull")
  log_k_values <- c(exponential = 0.0, weibull = 0.0)

  artifact_dir <- testthat::test_path("_artifacts", "mle_survival_gamma_compare")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE)
  }

  for (model_type in model_types) {
    log_k <- log_k_values[[model_type]]
    data <- if (model_type == "exponential") {
      mock_dataset_survival_exponential(
        n = n,
        var_v = var_v,
        var_epsilon = var_epsilon,
        beta_g = beta_g,
        beta_w = beta_w,
        theta = theta,
        censor_rate = censor_rate
      )
    } else {
      mock_dataset_survival_weibull(
        n = n,
        var_v = var_v,
        var_epsilon = var_epsilon,
        beta_g = beta_g,
        beta_w = beta_w,
        theta = theta,
        log_k = log_k,
        censor_rate = censor_rate
      )
    }

    start_beta <- rep(0, ncol(data$w) + 2)
    start_gamma <- rep(0, ncol(data$w) + 2)
    start_delta <- if (model_type == "weibull") c(log_k = 0) else numeric(0)

    old_time <- system.time({
      old_fit <- hapr_mle_survival(
        event_time = data$event_time,
        event_status = data$event_status,
        gc = data$gc,
        w = data$w,
        improvement_ratio = improvement_ratio,
        model_type = model_type,
        start_beta = start_beta,
        start_delta = start_delta,
        control = list(maxit = 150)
      )
    })

    new_time <- system.time({
      new_fit <- hapr_mle_survival_gamma(
        event_time = data$event_time,
        event_status = data$event_status,
        gc = data$gc,
        w = data$w,
        improvement_ratio = improvement_ratio,
        model_type = model_type,
        start_gamma = start_gamma,
        start_delta = start_delta,
        control = list(maxit = 150)
      )
    })

    runtime_old_ms <- old_time[["elapsed"]] * 1000
    runtime_new_ms <- new_time[["elapsed"]] * 1000

    beta_old <- old_fit$parameters$beta
    beta_new <- new_fit$parameters$beta

    aligned_names <- intersect(names(beta_old), names(beta_new))
    beta_old <- beta_old[aligned_names]
    beta_new <- beta_new[aligned_names]

    abs_diff <- abs(beta_old - beta_new)

    se_old <- if (is.null(old_fit$standard_errors)) {
      rep(NA_real_, length(aligned_names))
    } else {
      old_fit$standard_errors[aligned_names]
    }
    se_new <- if (is.null(new_fit$standard_errors)) {
      rep(NA_real_, length(aligned_names))
    } else {
      new_fit$standard_errors[aligned_names]
    }

    se_diff <- se_new - se_old

    comparison_table <- data.frame(
      Parameter = aligned_names,
      Beta_Old = beta_old,
      Beta_New = beta_new,
      Abs_Diff = abs_diff,
      SE_Old = se_old,
      SE_New = se_new,
      SE_Diff = se_diff,
      Runtime_Old_ms = rep(runtime_old_ms, length(aligned_names)),
      Runtime_New_ms = rep(runtime_new_ms, length(aligned_names)),
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    if (model_type == "weibull") {
      delta_old <- unlist(old_fit$parameters$delta)
      delta_new <- unlist(new_fit$parameters$delta)
      if (length(delta_old) == 1 && length(delta_new) == 1) {
        comparison_table <- rbind(
          comparison_table,
          data.frame(
            Parameter = names(delta_old),
            Beta_Old = delta_old,
            Beta_New = delta_new,
            Abs_Diff = abs(delta_old - delta_new),
            SE_Old = NA_real_,
            SE_New = NA_real_,
            SE_Diff = NA_real_,
            Runtime_Old_ms = runtime_old_ms,
            Runtime_New_ms = runtime_new_ms,
            row.names = NULL,
            stringsAsFactors = FALSE
          )
        )
      }
    }

    artifact_file <- file.path(artifact_dir, paste0("compare_", model_type, ".csv"))
    write.csv(comparison_table, artifact_file, row.names = FALSE)

    tol <- 1e-3
    expect_lt(max(abs_diff, na.rm = TRUE), tol)

    se_ok <- (se_new + 1e-8) >= se_old
    se_ok[is.na(se_ok)] <- TRUE
    expect_true(all(se_ok))

    expect_true(mean(se_diff, na.rm = TRUE) >= 0)
  }
})
