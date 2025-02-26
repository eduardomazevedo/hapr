regress_lm_y_on_gc_and_w <- function(y, gc, w) {
  # Combine gc with w into a data frame
  df <- data.frame(y = y, gc = gc, w)

  # Fit linear model using lm (for readability)
  model <- lm(y ~ ., data = df)
  model_summary <- summary(model)

  # Extract coefficients, standard errors, and variance-covariance matrix
  gamma <- coef(model)
  se_gamma <- model_summary$coefficients[, "Std. Error"]
  vcov_gamma <- vcov(model)

  # Compute R² for y ~ gc + w
  r2_gc_and_w <- model_summary$r.squared

  # Compute variance of error term
  var_error_y_on_gc_and_w <- model_summary$sigma^2

  # Fit regression model of y on gc only (including intercept)
  model_y_on_gc <- lm(y ~ gc, data = df)
  model_y_on_gc_summary <- summary(model_y_on_gc)
  r2_gc <- model_y_on_gc_summary$r.squared

  # Return results as a list
  return(list(
    gamma = gamma,
    se_gamma = se_gamma,
    vcov_gamma = vcov_gamma,
    r2_gc_and_w = r2_gc_and_w,
    r2_gc = r2_gc,
    var_error_y_on_gc_and_w = var_error_y_on_gc_and_w
  ))
}
