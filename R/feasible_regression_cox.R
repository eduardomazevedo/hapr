#' Regress survival outcome y on gc and w using Cox proportional hazards model
#' @param y Survival outcome (Surv object)
#' @param gc Polygenic risk score (normalized)
#' @param w Control variables data frame
#' @return A list containing:
#'   \item{gamma}{Regression coefficients}
#'   \item{se_gamma}{Standard errors of coefficients}
#'   \item{vcov_gamma}{Variance-covariance matrix}
#'   \item{cox_model}{Fitted Cox proportional hazards model}
#' @details
#' This function fits a Cox proportional hazards model of y on gc and w and returns the coefficients,
#' standard errors, variance-covariance matrix, and fitted model object.
#' Used in the first stage of hapr_cox.
#' @noRd
feasible_regression_cox <- function(y, gc, w) {
  # Combine gc with w into a data frame (y is passed separately)
  df <- data.frame(gc = gc, w)

  # Fit Cox proportional hazards model
  model <- survival::coxph(y ~ ., data = df)
  model_summary <- summary(model)

  # Extract coefficients, standard errors, and variance-covariance matrix
  gamma <- coef(model)
  se_gamma <- model_summary$coefficients[, "se(coef)"]
  vcov_gamma <- vcov(model)

  # Return results as a list
  list(
    gamma = gamma,
    se_gamma = se_gamma,
    vcov_gamma = vcov_gamma,
    cox_model = model
  )
}
