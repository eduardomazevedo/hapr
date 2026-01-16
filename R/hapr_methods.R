#' Print method for hapr_first_stage_fit objects
#'
#' Prints a concise summary of a hapr_first_stage_fit object.
#'
#' @param x A hapr_first_stage_fit object
#' @param ... Additional arguments (not used)
#'
#' @return The input object, invisibly
#' @export
print.hapr_first_stage_fit <- function(x, ...) {
  monkey <- "\U1F435"
  cat(monkey, " HAPR First Stage Fit ", monkey, "\n")
  cat("-------------------------------------------\n")
  
  cat("Model type:", x$model_type, "\n\n")
  
  cat("Theta coefficients (gc ~ w):\n")
  theta_to_show <- x$parameters$theta[1:min(5, length(x$parameters$theta))]
  coef_table <- data.frame(Estimate = theta_to_show, row.names = names(theta_to_show))
  print(coef_table, digits = 4)
  if (length(x$parameters$theta) > 5) {
    cat("  Showing first 5 of", length(x$parameters$theta), "coefficients\n")
  }
  cat("\n")
  
  cat("Gamma coefficients (y ~ gc + w):\n")
  gamma_to_show <- x$parameters$gamma[1:min(5, length(x$parameters$gamma))]
  coef_table <- data.frame(Estimate = gamma_to_show, row.names = names(gamma_to_show))
  print(coef_table, digits = 4)
  if (length(x$parameters$gamma) > 5) {
    cat("  Showing first 5 of", length(x$parameters$gamma), "coefficients\n")
  }
  cat("\n")
  
  cat("Var(v + epsilon):", sprintf("%.4f", x$parameters$var_v_plus_var_epsilon), "\n")
  cat("Max improvement ratio:", sprintf("%.4f", x$stats$max_improvement_ratio), "\n")
  cat("Var(w*theta):", sprintf("%.4f", x$stats$var_wtheta), "\n")
  
  invisible(x)
}

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
  beta_to_show <- x$parameters$beta[1:min(5, length(x$parameters$beta))]
  coef_table <- data.frame(Estimate = beta_to_show, row.names = names(beta_to_show))
  print(coef_table, digits = 4)
  if (length(x$parameters$beta) > 5) {
    cat("  Showing first 5 of", length(x$parameters$beta), "coefficients\n")
  }
  cat("\n")

  cat("Improvement ratio:", sprintf("%.4f", x$stats$improvement_ratio), "\n")

  if (!is.na(x$stats$r2_current)) {
    cat("R\U00B2 current:", sprintf("%.4f", x$stats$r2_current), "\n")
  }

  if (!is.na(x$stats$r2_future)) {
    cat("R\U00B2 future:", sprintf("%.4f", x$stats$r2_future), "\n")
  }

  cat("Max improvement ratio:", sprintf("%.4f", x$stats$max_improvement_ratio), "\n")

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
#' @export
summary.hapr_fit <- function(object, ...) {
  result <- list(
    model_type = object$model_type,
    beta = object$parameters$beta,
    sd_beta = object$standard_errors,
    ci_beta = object$ci_beta,
    gamma = object$parameters$gamma,
    theta = object$parameters$theta,
    var_v = object$stats$var_v,
    var_epsilon = object$stats$var_epsilon,
    improvement_ratio = object$stats$improvement_ratio,
    max_improvement_ratio = object$stats$max_improvement_ratio,
    r2_current = object$stats$r2_current,
    r2_future = object$stats$r2_future,
    posterior = object$stats$posterior
  )
  
  # Add model-specific parameters if present
  if (!is.null(object$parameters$var_eta)) {
    result$var_eta <- object$parameters$var_eta
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
  print(data.frame(Estimate = x$beta, row.names = names(x$beta)), digits = 4)
  cat("\n")

  if (!is.null(x$sd_beta)) {
    cat("Standard errors (delta method):\n")
    print(data.frame(Std.Error = x$sd_beta, row.names = names(x$sd_beta)), digits = 4)
    cat("\n")
  }

  if (!is.null(x$ci_beta)) {
    cat("95% Confidence Intervals for Beta (delta method):\n")
    print(x$ci_beta, digits = 4)
    cat("\n")
  }

  cat("Gamma coefficients (current PRS effects):\n")
  print(data.frame(Estimate = x$gamma, row.names = names(x$gamma)), digits = 4)
  cat("\n")

  cat("Theta coefficients (PRS decomposition):\n")
  print(data.frame(Estimate = x$theta, row.names = names(x$theta)), digits = 4)
  cat("\n")

  cat("Model Statistics:\n")
  cat("-----------------\n")
  cat("Improvement ratio:", sprintf("%.4f", x$improvement_ratio), "\n")
  cat("Max improvement ratio:", sprintf("%.4f", x$max_improvement_ratio), "\n")
  if (!is.na(x$r2_current)) {
    cat("R\U00B2 current:", sprintf("%.4f", x$r2_current), "\n")
  }
  if (!is.na(x$r2_future)) {
    cat("R\U00B2 future:", sprintf("%.4f", x$r2_future), "\n")
  }

  cat("Var(v):", sprintf("%.4f", x$var_v), "\n")
  cat("Var(epsilon):", sprintf("%.4f", x$var_epsilon), "\n\n")

  cat("Posterior Parameters:\n")
  cat("--------------------\n")
  cat("a:", sprintf("%.4f", x$posterior$a), "\n")
  cat("b:", sprintf("%.4f", x$posterior$b), "\n")
  cat("c:", sprintf("%.4f", x$posterior$c), "\n\n")

  if (!is.null(x$var_eta)) {
    cat("Var(eta):", sprintf("%.4f", x$var_eta), "\n\n")
  }
  
  invisible(x)
}
