#' Regress binary y on gc and w using probit regression
#' @param y Binary outcome variable (factor)
#' @param gc Polygenic risk score (normalized)
#' @param w Control variables data frame
#' @return A list containing:
#'   \item{gamma}{Regression coefficients}
#'   \item{se_gamma}{Standard errors of coefficients}
#'   \item{vcov_gamma}{Variance-covariance matrix}
#'   \item{r2_gc_and_w}{Liability R² for y ~ gc + w}
#'   \item{r2_gc}{Liability R² for y ~ gc}
#'   \item{var_error_y_on_gc_and_w}{Residual variance}
#' @details
#' This function fits a probit regression model of y on gc and w and returns the coefficients,
#' standard errors, variance-covariance matrix, liability R² values, and residual variance.
#' Used in the first stage of hapr_probit.
#' @noRd
feasible_regression_probit <- function(y, gc, w) {
  # Combine gc with w into a data frame
  df <- data.frame(y = y, gc = gc, w)

  # Fit linear model using lm (for readability)
  model <- glm(y ~ ., data = df, family = binomial(link = "probit"))

  # Extract coefficients, standard errors, and variance-covariance matrix
  gamma <- coef(model)
  se_gamma <- summary(model)$coefficients[, "Std. Error"]
  vcov_gamma <- vcov(model)

  # Compute R² for y ~ gc + w
  r2_gc_and_w <- r2_liability_probit(model)

  # Compute variance of error term
  var_error_y_on_gc_and_w <- summary(model)$sigma^2

  # Fit regression model of y on gc only (including intercept)
  model_y_on_gc <- glm(y ~ gc, data = df, family = binomial(link = "probit"))
  r2_gc <- r2_liability_probit(model_y_on_gc)

  # Return results as a list
  list(
    gamma = gamma,
    se_gamma = se_gamma,
    vcov_gamma = vcov_gamma,
    r2_gc_and_w = r2_gc_and_w,
    r2_gc = r2_gc,
    var_error_y_on_gc_and_w = var_error_y_on_gc_and_w
  )
}
