#' Print method for hapr_fit objects
#'
#' Prints a concise summary of a hapr_fit object.
#'
#' @param x A hapr_fit object
#' @param ... Additional arguments (not used)
#'
#' @return The input object, invisibly
#' @export
print.hapr_fit <- function(x, ...) {
  monkey <- "\U1F435"
  cat(monkey, " HAPR (Heritability Adjusted Prediction) Model ", monkey, "\n")
  cat("-------------------------------------------\n")

  cat("Model type:", x$model_type, "\n\n")

  cat("Beta coefficients (future PRS effects):\n")
  beta_to_show <- x$coefficients$beta[1:min(5, length(x$coefficients$beta))]
  coef_table <- data.frame(Estimate = beta_to_show)
  print(coef_table, digits = 4)
  if (length(x$coefficients$beta) > 5) {
    cat("  Showing first 5 of", length(x$coefficients$beta), "coefficients\n")
  }
  cat("\n")

  cat("Improvement ratio:", sprintf("%.4f", x$stats$improvement_ratio),
      "(", x$stats$heritability_source, ")\n")

  if (!is.na(x$stats$r2_current)) {
    cat("R² current:", sprintf("%.4f", x$stats$r2_current),
        "(", x$stats$r2_current_source, ")\n")
  }

  if (!is.na(x$stats$r2_future)) {
    cat("R² future:", sprintf("%.4f", x$stats$r2_future), "\n")
  }

  cat("Max improvement ratio:", sprintf("%.4f", x$stats$max_improvement_ratio), "\n")

  if (x$model_type == "cox") {
    if (!is.null(x$additional_parameters$base_hazard_conversion_ratio)) {
      cat("Base hazard conversion ratio:",
          sprintf("%.4f", x$additional_parameters$base_hazard_conversion_ratio), "\n")
    }
  }

  invisible(x)
}

#' Summary method for hapr_fit objects
#'
#' Provides a detailed summary of a hapr_fit object.
#'
#' @param object A hapr_fit object
#' @param ... Additional arguments (not used)
#'
#' @return A summary.hapr_fit object containing the summary information
#'
#' @details
#' For Cox models, the reported standard errors and confidence intervals are based on delta method approximations.
#' These may underestimate true uncertainty, especially for the `gf` term, due to:
#' \itemize{
#'   \item The use of a partial likelihood (not full likelihood)
#'   \item Nonlinear transformation of parameters not fully captured by the Jacobian
#' }
#' For production inference, consider resampling methods such as bootstrap.
#' @export
summary.hapr_fit <- function(object, ...) {
  result <- list(
    model_type = object$model_type,
    beta = object$coefficients$beta,
    sd_beta = object$standard_errors,
    ci_beta = object$ci_beta,
    gamma = object$coefficients$gamma,
    theta = object$coefficients$theta,
    var_v = object$stats$var_v,
    var_epsilon = object$stats$var_epsilon,
    improvement_ratio = object$stats$improvement_ratio,
    max_improvement_ratio = object$stats$max_improvement_ratio,
    r2_current = object$stats$r2_current,
    r2_future = object$stats$r2_future,
    heritability_source = object$stats$heritability_source,
    r2_current_source = object$stats$r2_current_source,
    posterior = object$stats$posterior
  )

  if (object$model_type == "cox") {
    result$base_hazard_conversion_ratio <- object$additional_parameters$base_hazard_conversion_ratio
  }

  class(result) <- "summary.hapr_fit"
  result
}

#' @export
print.summary.hapr_fit <- function(x, ...) {
  monkey <- "\U1F435"
  cat(monkey, " HAPR (Heritability Adjusted Prediction) Model Summary ", monkey, "\n")
  cat("=================================================\n\n")

  cat("Model type:", x$model_type, "\n\n")

  cat("Beta coefficients (future PRS effects):\n")
  print(data.frame(Estimate = x$beta), digits = 4)
  cat("\n")

  if (!is.null(x$sd_beta)) {
    cat("Standard errors (delta method):\n")
    print(data.frame(Std.Error = x$sd_beta), digits = 4)
    cat("\n")
  }

  if (!is.null(x$ci_beta)) {
    cat("95% Confidence Intervals for Beta (delta method):\n")
    print(x$ci_beta, digits = 4)
    cat("\n")
  }

  cat("Gamma coefficients (current PRS effects):\n")
  print(data.frame(Estimate = x$gamma), digits = 4)
  cat("\n")

  cat("Theta coefficients (PRS decomposition):\n")
  print(data.frame(Estimate = x$theta), digits = 4)
  cat("\n")

  cat("Model Statistics:\n")
  cat("-----------------\n")
  cat("Improvement ratio:", sprintf("%.4f", x$improvement_ratio),
      "(", x$heritability_source, ")\n")
  cat("Max improvement ratio:", sprintf("%.4f", x$max_improvement_ratio), "\n")
  if (!is.na(x$r2_current)) {
    cat("R² current:", sprintf("%.4f", x$r2_current),
        "(", x$r2_current_source, ")\n")
  }
  if (!is.na(x$r2_future)) {
    cat("R² future:", sprintf("%.4f", x$r2_future), "\n")
  }

  cat("Var(v):", sprintf("%.4f", x$var_v), "\n")
  cat("Var(epsilon):", sprintf("%.4f", x$var_epsilon), "\n\n")

  cat("Posterior Parameters:\n")
  cat("--------------------\n")
  cat("a:", sprintf("%.4f", x$posterior$a), "\n")
  cat("b:", sprintf("%.4f", x$posterior$b), "\n")
  cat("c:", sprintf("%.4f", x$posterior$c), "\n\n")

  if (x$model_type == "cox") {
    cat("Cox Model-Specific Statistics:\n")
    cat("----------------------------\n")
    if (!is.null(x$base_hazard_conversion_ratio)) {
      cat("Base Hazard Conversion Ratio:", sprintf("%.4f", x$base_hazard_conversion_ratio), "\n")
      cat("Baseline hazard: [matrix omitted for brevity]\n")
    }
  }
  invisible(x)
}
