#' Extract baseline hazard from a HAPR Cox model
#'
#' @param fit A hapr_fit object from a Cox model
#' @param covariates Character string indicating which covariates to use, one of: 'gf_w' (default), 'gc_w', or 'w'
#' @return A data frame containing the baseline hazard with time and hazard columns
#' @export
hapr_basehaz <- function(fit, covariates = "gf_w") {
  if (!inherits(fit, "hapr_fit")) stop("hapr_basehaz function only works for hapr_fit objects.")
  if (fit$model_type != "cox") stop("hapr_basehaz function only works for Cox models.")
  if (!(covariates %in% c("gf_w", "gc_w", "w"))) stop("covariates must be 'gf_w', 'gc_w', or 'w'.")
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
hapr_survfit <- function(fit,
                         covariates = "gf_w",
                         newdata,
                         start.time = NULL,
                         aggregate = NULL,
                         conf.int = FALSE,
                         conf.level = 0.95,
                         conf.method = c("bootstrap"),
                         n.boot = 200) {
  if (!inherits(fit, "hapr_fit")) stop("hapr_survfit requires a hapr_fit object.")
  if (fit$model_type != "cox") stop("hapr_survfit only works for Cox models.")
  if (!(covariates %in% c("gf_w", "gc_w", "w"))) stop("covariates must be 'gf_w', 'gc_w', or 'w'.")

  conf.method <- match.arg(conf.method)
  baseline_hazard <- hapr_basehaz(fit, covariates)
  relative_risk <- predict.hapr_fit(fit, newdata, covariates, type = "risk")[[paste0("y_hat_", covariates)]]

  cumhaz <- outer(baseline_hazard$hazard, relative_risk, "*")
  surv <- exp(-cumhaz)
  time <- baseline_hazard$time

  if (!is.null(start.time)) {
    if (!is.numeric(start.time) || length(start.time) != 1) stop("start.time must be a single numeric value.")
    keep_idx <- which(time >= start.time)
    if (length(keep_idx) == 0) stop("start.time is beyond the range of survival times.")
    t0_idx <- max(which(time < start.time))
    haz0 <- baseline_hazard$hazard[t0_idx]
    surv_t0 <- exp(-haz0 * relative_risk)
    time <- time[keep_idx]
    cumhaz <- cumhaz[keep_idx, , drop = FALSE]
    surv <- surv[keep_idx, , drop = FALSE]
    surv <- surv / matrix(surv_t0, nrow = length(time), ncol = length(surv_t0), byrow = TRUE)
  }

  result <- list()
  result$time <- time
  result$relative_risk <- relative_risk
  result$surv <- surv

  if (!is.null(aggregate)) {
    if (is.expression(substitute(aggregate))) {
      aggregate <- eval(substitute(aggregate), envir = newdata)
    }
    if (!is.logical(aggregate) || length(aggregate) != nrow(newdata)) {
      stop("aggregate must be a logical vector matching rows in newdata.")
    }
    idx <- which(aggregate)
    if (length(idx) == 0) stop("No subjects match the aggregation condition.")
    surv_subset <- surv[, idx, drop = FALSE]
    surv_avg <- rowMeans(surv_subset)
    result$surv_avg <- surv_avg
    result$aggregate_idx <- idx

    if (conf.int && conf.method == "bootstrap") {
      boot_mat <- matrix(NA, nrow = length(time), ncol = n.boot)
      set.seed(123)
      for (b in seq_len(n.boot)) {
        samp <- sample(idx, replace = TRUE)
        boot_mat[, b] <- rowMeans(surv[, samp, drop = FALSE])
      }
      alpha <- 1 - conf.level
      result$surv_avg_lower <- apply(boot_mat, 1, quantile, probs = alpha / 2)
      result$surv_avg_upper <- apply(boot_mat, 1, quantile, probs = 1 - alpha / 2)
    }
  }

  class(result) <- "hapr_survfit"
  return(result)
}

#' Plot survival curves from a HAPR survival fit
#'
#' @param x A hapr_survfit object returned by survfit.hapr_fit
#' @param mode Character string indicating how to select curves for plotting:
#'        "percentiles" (default) to show curves at specific risk percentiles
#'        "subjects" to show curves for each individual subject
#' @param percentiles Numeric vector of probability values (0-1) indicating which
#'        percentiles to plot when mode="percentiles"
#' @param ... Additional arguments (not currently used)
#' @return A plot of the survival curves (using ggplot2 if available, otherwise base R)
#' @export
plot.hapr_survfit <- function(x,
                              mode = "percentiles",
                              percentiles = c(0.01, 0.10, 0.22, 0.50, 0.75, 0.90, 0.99),
                              ...) {
  time <- x$time
  surv_matrix <- x$surv
  relative_risk <- x$relative_risk
  num_subjects <- length(relative_risk)
  if (num_subjects == 0) stop("No subjects available for plotting.")
  plot_data <- data.frame(time = numeric(), survival = numeric(), label = character())

  if (mode == "percentiles") {
    quantile_thresholds <- stats::quantile(relative_risk, probs = percentiles, na.rm = TRUE)
    for (p in seq_along(percentiles)) {
      percentile_value <- quantile_thresholds[p]
      closest_idx <- which.min(abs(relative_risk - percentile_value))
      if (length(closest_idx) == 0 || closest_idx > num_subjects) next
      surv_curve <- surv_matrix[, closest_idx]
      temp_data <- data.frame(time = time, survival = surv_curve,
                              label = paste0("Percentile: ", percentiles[p] * 100, "%"))
      plot_data <- rbind(plot_data, temp_data)
    }
  } else if (mode == "subjects") {
    for (i in seq_along(relative_risk)) {
      temp_data <- data.frame(time = time, survival = surv_matrix[, i], label = paste0("Subject ", i))
      plot_data <- rbind(plot_data, temp_data)
    }
  } else stop("Invalid mode. Use 'percentiles' or 'subjects'.")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = survival, color = label)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::labs(title = "Survival Curves", x = "Time", y = "Survival Probability") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::scale_color_viridis_d() +
      ggplot2::theme(aspect.ratio = 9 / 16)
  } else {
    unique_labels <- unique(plot_data$label)
    colors <- rainbow(length(unique_labels))
    plot(NULL, xlim = range(time), ylim = c(0, 1), xlab = "Time", ylab = "Survival Probability",
         main = "Survival Curves")
    for (i in seq_along(unique_labels)) {
      subset_data <- subset(plot_data, label == unique_labels[i])
      lines(subset_data$time, subset_data$survival, col = colors[i], lwd = 2)
    }
    legend("topright", legend = unique_labels, col = colors, lwd = 2, bty = "n")
  }
}
