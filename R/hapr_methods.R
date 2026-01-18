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
  
  print_ci_table("Theta coefficients (gc ~ w):",
                 x$parameters$theta,
                 x$vcov_parameters$theta)
  
  if (!is.null(x$parameters$gamma) && !is.null(x$vcov_parameters$gamma)) {
    print_ci_table("Gamma coefficients (y ~ gc + w):",
                   x$parameters$gamma,
                   x$vcov_parameters$gamma)
  } else {
    cat("Gamma coefficients (y ~ gc + w):\n")
    cat("  (not estimated for model_type = 'mle')\n\n")
  }
  
  print_ci_table("Var(v + epsilon):",
                 c(var_v_plus_var_epsilon = x$parameters$var_v_plus_var_epsilon),
                 x$vcov_parameters$var_v_plus_var_epsilon)
  cat("Max improvement ratio:", sprintf("%.4f", x$stats$max_improvement_ratio), "\n")
  cat("Var(w*theta):", sprintf("%.4f", x$stats$var_wtheta), "\n")
  
  invisible(x)
}

print_coef_table <- function(values, max_show = 5) {
  n_coef <- nrow(values)
  n_to_show <- min(max_show, n_coef)
  table <- values[1:n_to_show, , drop = FALSE]
  list(table = table, truncated = n_coef > max_show, total = n_coef)
}

print_named_values <- function(title, values, digits = 4) {
  cat(title, "\n")
  table <- data.frame(Estimate = values, row.names = names(values))
  print(table, digits = digits)
  cat("\n")
}

make_ci_table <- function(estimates, se, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(
    Estimate = estimates,
    Std.Error = se,
    Lower = estimates - z * se,
    Upper = estimates + z * se,
    row.names = names(estimates),
    check.names = FALSE
  )
}

delta_list_to_vector <- function(delta) {
  if (is.null(delta)) {
    return(numeric(0))
  }
  if (is.list(delta)) {
    return(unlist(delta, use.names = TRUE))
  }
  delta
}

extract_se <- function(estimates, vcov) {
  if (is.null(vcov)) {
    return(NULL)
  }
  if (is.matrix(vcov)) {
    se <- sqrt(diag(vcov))
    names(se) <- names(estimates)
    return(se)
  }
  se <- sqrt(as.numeric(vcov))
  names(se) <- names(estimates)
  se
}

subset_vcov <- function(vcov, idx) {
  if (is.null(vcov) || !is.matrix(vcov) || length(idx) == 0) {
    return(NULL)
  }
  vcov[idx, idx, drop = FALSE]
}

print_ci_table <- function(title, estimates, vcov, max_show = 5) {
  se <- extract_se(estimates, vcov)
  if (is.null(se)) {
    print_named_values(title, estimates)
    return(invisible(NULL))
  }
  cat(title, "\n")
  ci <- make_ci_table(estimates, se)
  coef_table <- print_coef_table(ci, max_show = max_show)
  print(coef_table$table, digits = 4)
  if (coef_table$truncated) {
    cat("  Showing first 5 of", coef_table$total, "coefficients\n")
  }
  cat("\n")
  invisible(ci)
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

  print_ci_table("Beta coefficients (future PRS effects):",
                 x$parameters$beta,
                 x$vcov_parameters$beta)
  cat("Note: Standard errors are delta-method approximations and may be conservative.\n\n")

  if (!is.null(x$parameters$theta) && !is.null(x$vcov_parameters$theta)) {
    print_ci_table("Theta coefficients (gc ~ w):",
                   x$parameters$theta,
                   x$vcov_parameters$theta)
  }

  if (!is.null(x$regressions$y_on_gc_w$sigma_squared) &&
      !is.null(x$regressions$y_on_gc_w$var_sigma_squared)) {
    print_ci_table("Delta parameter (sigma^2_y):",
                   c(sigma2_y = x$regressions$y_on_gc_w$sigma_squared),
                   x$regressions$y_on_gc_w$var_sigma_squared)
  }

  if (!is.null(x$regressions$gc_on_w$sigma_squared) &&
      !is.null(x$regressions$gc_on_w$var_sigma_squared)) {
    print_ci_table("Stage 1 variance (v + epsilon):",
                   c(var_v_plus_var_epsilon = x$regressions$gc_on_w$sigma_squared),
                   x$regressions$gc_on_w$var_sigma_squared)
  }

  cat("Improvement ratio:", sprintf("%.4f", x$stats$improvement_ratio), "\n")

  if (!is.na(x$stats$r2_current)) {
    # Use scientific notation for very small values (< 0.0001)
    if (abs(x$stats$r2_current) < 0.0001 && x$stats$r2_current != 0) {
      cat("R\U00B2 current:", sprintf("%.4e", x$stats$r2_current), "\n")
    } else {
      cat("R\U00B2 current:", sprintf("%.4f", x$stats$r2_current), "\n")
    }
  }

  if (!is.na(x$stats$r2_future)) {
    # Use scientific notation for very small values (< 0.0001)
    if (abs(x$stats$r2_future) < 0.0001 && x$stats$r2_future != 0) {
      cat("R\U00B2 future:", sprintf("%.4e", x$stats$r2_future), "\n")
    } else {
      cat("R\U00B2 future:", sprintf("%.4f", x$stats$r2_future), "\n")
    }
  }

  cat("Max improvement ratio:", sprintf("%.4f", x$stats$max_improvement_ratio), "\n")

  invisible(x)
}

#' Print method for hapr_mle_fit objects
#'
#' Prints a concise summary of a hapr_mle_fit object.
#'
#' @param x A hapr_mle_fit object
#' @param ... Additional arguments (not used)
#'
#' @return The input object, invisibly
#' @export
print.hapr_mle_fit <- function(x, ...) {
  monkey <- "\U1F435"
  cat(monkey, " HAPR MLE Fit ", monkey, "\n")
  cat("-------------------------------------------\n")

  cat("Model type:", x$model_type, "\n\n")

  order_beta <- x$vcov_parameters$order$beta
  if (is.null(order_beta)) {
    order_beta <- seq_along(x$parameters$beta)
  }
  vcov_beta <- subset_vcov(x$vcov_parameters$all, order_beta)
  if (!is.null(vcov_beta)) {
    print_ci_table("Beta coefficients:", x$parameters$beta, vcov_beta)
    cat("Note: MLE standard errors ignore first-stage uncertainty.\n\n")
  } else {
    print_named_values("Beta coefficients:",
                       x$parameters$beta[1:min(5, length(x$parameters$beta))])
    if (length(x$parameters$beta) > 5) {
      cat("  Showing first 5 of", length(x$parameters$beta), "coefficients\n\n")
    }
  }

  delta_values <- delta_list_to_vector(x$parameters$delta)
  if (length(delta_values) > 0) {
    order_delta <- x$vcov_parameters$order$delta
    if (is.null(order_delta)) {
      order_delta <- seq_along(delta_values)
    }
    vcov_delta <- subset_vcov(x$vcov_parameters$all, order_delta)
    if (!is.null(vcov_delta)) {
      print_ci_table("Delta parameters:", delta_values, vcov_delta)
    } else {
      print_named_values("Delta parameters:", delta_values)
    }
  }

  if (!is.null(x$parameters$var_v_plus_var_epsilon) &&
      !is.null(x$regressions$gc_on_w$var_sigma_squared)) {
    print_ci_table("Stage 1 variance (v + epsilon):",
                   c(var_v_plus_var_epsilon = x$parameters$var_v_plus_var_epsilon),
                   x$regressions$gc_on_w$var_sigma_squared)
  }

  if (!is.null(x$parameters$theta) &&
      !is.null(x$regressions$gc_on_w$vcov_coefficients)) {
    print_ci_table("Theta coefficients (gc ~ w):",
                   x$parameters$theta,
                   x$regressions$gc_on_w$vcov_coefficients)
  }

  cat("Improvement ratio:", sprintf("%.4f", x$stats$improvement_ratio), "\n")
  cat("Max improvement ratio:", sprintf("%.4f", x$stats$max_improvement_ratio), "\n")
  cat("Var(v):", sprintf("%.4f", x$stats$var_v), "\n")
  cat("Var(epsilon):", sprintf("%.4f", x$stats$var_epsilon), "\n")

  if (!is.null(x$opt)) {
    cat("Convergence:", x$opt$convergence, "\n")
    if (!is.null(x$opt$value)) {
      cat("Neg. log-likelihood:", sprintf("%.4f", x$opt$value), "\n")
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
    cat("Note: Standard errors are delta-method approximations and may be conservative.\n\n")
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
    # Use scientific notation for very small values (< 0.0001)
    if (abs(x$r2_current) < 0.0001 && x$r2_current != 0) {
      cat("R\U00B2 current:", sprintf("%.4e", x$r2_current), "\n")
    } else {
      cat("R\U00B2 current:", sprintf("%.4f", x$r2_current), "\n")
    }
  }
  if (!is.na(x$r2_future)) {
    # Use scientific notation for very small values (< 0.0001)
    if (abs(x$r2_future) < 0.0001 && x$r2_future != 0) {
      cat("R\U00B2 future:", sprintf("%.4e", x$r2_future), "\n")
    } else {
      cat("R\U00B2 future:", sprintf("%.4f", x$r2_future), "\n")
    }
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
