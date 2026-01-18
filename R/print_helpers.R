#' Print helpers for HAPR objects
#'
#' @noRd
NULL

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

delta_list_to_vector <- function(delta) {
  if (is.null(delta)) {
    return(numeric(0))
  }
  if (is.list(delta)) {
    return(unlist(delta, use.names = TRUE))
  }
  delta
}
