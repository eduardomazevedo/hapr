#' Extract baseline hazard from a HAPR Cox model
#'
#' @param fit A hapr_fit object from a Cox model
#' @param covariates Character string indicating which covariates to use, one of: 'gf_w' (default), 'gc_w', or 'w'
#' @return A data frame containing the baseline hazard with time and hazard columns
#' @export
basehaz.hapr_fit <- function(fit, covariates = 'gf_w') {
  if (fit$model_type != "cox") {
    stop("Basehaz function only works for Cox models.")
  }
  if (!(covariates %in% c('gf_w', 'gc_w', 'w'))) {
    stop("Basehaz function only works for gf_w, gc_w or w.")
  }

  model_name <- paste0("y_on_", covariates)
  if (!model_name %in% names(fit$regressions)) {
    stop(sprintf("No regression found for covariate '%s'; skipping.", covariates))
  }

  fit$regressions[[model_name]]$baseline_hazard
}

#' Create a survival fit from a HAPR Cox model
#'
#' @param fit A hapr_fit object from a Cox model
#' @param covariates Character string indicating which covariates to use, one of: 'gf_w' (default), 'gc_w', or 'w'
#' @param newdata A data frame containing new observations for which to compute survival curves
#' @return A hapr_survfit object containing survival probabilities over time
#' @export
survfit.hapr_fit <- function(fit, covariates = 'gf_w', newdata) {
  if (fit$model_type != "cox") {
    stop("Basehaz function only works for Cox models.")
  }
  if (!(covariates %in% c('gf_w', 'gc_w', 'w'))) {
    stop("Basehaz function only works for gf_w, gc_w or w.")
  }

  baseline_hazard <- basehaz.hapr_fit(fit, covariates)
  relative_risk <- predict.hapr_fit(fit, newdata, covariates, type = "risk")
  column_name <- paste0("y_hat_", covariates)
  relative_risk <- relative_risk[[column_name]]

  result <- list()
  result$time <- baseline_hazard$time
  result$relative_risk <- relative_risk
  result$cumhaz <- outer(baseline_hazard$hazard, relative_risk, "*")
  result$surv <- exp(-result$cumhaz)

  class(result) <- "hapr_survfit"

  return(result)
}

#' Plot survival curves from a HAPR survival fit
#'
#' @param hapr_survfit_object A hapr_survfit object returned by survfit.hapr_fit
#' @param mode Character string indicating how to select curves for plotting:
#'        "percentiles" (default) to show curves at specific risk percentiles
#'        "subjects" to show curves for each individual subject
#' @param percentiles Numeric vector of probability values (0-1) indicating which
#'        percentiles to plot when mode="percentiles"
#' @return A ggplot2 object showing the survival curves
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal theme element_blank scale_color_viridis_d
#' @export
plot.hapr_survfit <- function(hapr_survfit_object,
                              mode = "percentiles",
                              percentiles = c(0.01, 0.10, 0.22, 0.50, 0.75, 0.90, 0.99)) {
  # Extract relevant data
  time <- hapr_survfit_object$time
  surv_matrix <- hapr_survfit_object$surv
  relative_risk <- hapr_survfit_object$relative_risk

  # Ensure survival matrix has subjects (rows: time, columns: subjects)
  num_subjects <- length(relative_risk)
  if (num_subjects == 0) {
    stop("No subjects available for plotting.")
  }

  plot_data <- data.frame(time = numeric(), survival = numeric(), label = character())

  if (mode == "percentiles") {
    # Compute percentile thresholds
    quantile_thresholds <- stats::quantile(relative_risk, probs = percentiles, na.rm = TRUE, type = 7)

    for (p in seq_along(percentiles)) {
      percentile_value <- quantile_thresholds[p]

      # Find the subject whose relative risk is closest to this percentile
      closest_idx <- which.min(abs(relative_risk - percentile_value))

      if (length(closest_idx) == 0 || closest_idx > num_subjects) {
        warning(sprintf("No valid subject found for percentile %s.", percentiles[p] * 100))
        next
      }

      surv_curve <- surv_matrix[, closest_idx]

      temp_data <- data.frame(
        time = time,
        survival = surv_curve,
        label = paste0("Percentile: ", percentiles[p] * 100, "%")
      )

      # Append new rows to plot_data
      plot_data <- rbind(plot_data, temp_data)
    }

    # Ensure we have data to plot
    if (nrow(plot_data) == 0) {
      stop("No valid survival curves found for the selected percentiles.")
    }
  } else if (mode == "subjects") {
    for (i in seq_along(relative_risk)) {
      temp_data <- data.frame(
        time = time,
        survival = surv_matrix[, i],
        label = paste0("Subject ", i)
      )

      # Append new rows to plot_data
      plot_data <- rbind(plot_data, temp_data)
    }
  } else {
    stop("Invalid mode. Use 'percentiles' or 'subjects'.")
  }

  # Create the survival plot with a 16:9 aspect ratio
  ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = survival, color = label)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::labs(title = "Survival Curves", x = "Time", y = "Survival Probability") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::scale_color_viridis_d() +
    ggplot2::theme(
      aspect.ratio = 9 / 16  # Set 16:9 aspect ratio
    )
}
