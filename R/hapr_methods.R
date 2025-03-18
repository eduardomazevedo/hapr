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
  monkey <- "\U1F435" # Happy monkey emoji
  
  cat(monkey, " HAPR (Heritability Adjusted Prediction) Model ", monkey, "\n")
  cat("-------------------------------------------\n")
  
  cat("Model type:", x$model_type, "\n\n")
  
  # Print beta coefficients (first 5 at most)
  cat("Beta coefficients (future PRS effects):\n")
  beta_to_show <- x$coefficients$beta[1:min(5, length(x$coefficients$beta))]
  coef_table <- data.frame(Estimate = beta_to_show)
  print(coef_table, digits = 4)
  if (length(x$coefficients$beta) > 5) {
    cat("  Showing first 5 of", length(x$coefficients$beta), "coefficients\n")
  }
  cat("\n")
  
  # Key statistics
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
  
  # Cox-specific information without printing baseline hazard
  if (x$model_type == "cox" && !is.null(x$additional_parameters$base_hazard_conversion_ratio)) {
    cat("Base hazard conversion ratio:", 
        sprintf("%.4f", x$additional_parameters$base_hazard_conversion_ratio), "\n")
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
#' @export
summary.hapr_fit <- function(object, ...) {
  result <- list(
    model_type = object$model_type,
    beta = object$coefficients$beta,
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
  
  # Add model-specific information
  if (object$model_type == "cox" && !is.null(object$additional_parameters$base_hazard_conversion_ratio)) {
    result$base_hazard_conversion_ratio = object$additional_parameters$base_hazard_conversion_ratio
  }
  
  class(result) <- "summary.hapr_fit"
  result
}

#' Print method for summary.hapr_fit objects
#'
#' @param x A summary.hapr_fit object
#' @param ... Additional arguments (not used)
#'
#' @return The input object, invisibly
#' @export
print.summary.hapr_fit <- function(x, ...) {
  monkey <- "\U1F435" # Happy monkey emoji
  
  cat(monkey, " HAPR (Heritability Adjusted Prediction) Model Summary ", monkey, "\n")
  cat("=================================================\n\n")
  
  cat("Model type:", x$model_type, "\n\n")
  
  # Print coefficients in the specified order
  cat("Beta coefficients (future PRS effects):\n")
  print(data.frame(Estimate = x$beta), digits = 4)
  cat("\n")
  
  cat("Gamma coefficients (current PRS effects):\n")
  print(data.frame(Estimate = x$gamma), digits = 4)
  cat("\n")
  
  cat("Theta coefficients (PRS decomposition):\n")
  print(data.frame(Estimate = x$theta), digits = 4)
  cat("\n")
  
  # Statistics
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
  
  # Posterior parameters
  cat("Posterior Parameters:\n")
  cat("--------------------\n")
  cat("a:", sprintf("%.4f", x$posterior$a), "\n")
  cat("b:", sprintf("%.4f", x$posterior$b), "\n")
  cat("c:", sprintf("%.4f", x$posterior$c), "\n\n")
  
  # Cox-specific information
  if (x$model_type == "cox" && !is.null(x$base_hazard_conversion_ratio)) {
    cat("Cox Model-Specific Statistics:\n")
    cat("----------------------------\n")
    cat("Base Hazard Conversion Ratio:", 
        sprintf("%.4f", x$base_hazard_conversion_ratio), "\n")
    cat("Baseline hazard: [matrix omitted for brevity]\n\n")
  }
  
  invisible(x)
} 