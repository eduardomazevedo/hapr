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
