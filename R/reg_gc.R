regress_gc_on_w <- function(gc, w) {
  # Convert inputs into a data frame for readability
  df <- data.frame(gc = gc, w)

  # Fit the linear model using lm (for readability)
  model_gc <- lm(gc ~ ., data = df)
  model_gc_summary <- summary(model_gc)
  # Extract coefficients, standard errors, and variance-covariance matrix
  theta <- coef(model_gc)
  se_theta <- model_gc_summary$coefficients[, "Std. Error"]
  vcov_theta <- vcov(model_gc)

  # Extract unexplained variance (residual variance)
  var_residual <- sigma(model_gc)^2
  var_total <- var_residual # This is var_v + var_epsilon that gets saved

  # Compute explained variance (R² equivalent)
  var_wtheta <- var(fitted(model_gc))

  # Dearling with numerical issues
  # Normalize variances to sum to 1 if needed
  var_sum <- var_residual + var_wtheta
  if (var_sum > 0) {
    var_total <- var_residual / var_sum
    var_wtheta <- var_wtheta / var_sum
  } else {
    warning("Total variance sum is zero or negative. Check inputs.")
    var_total <- NA
    var_wtheta <- NA
  }


  # Compute maximum improvement ratio
  max_improvement_ratio <- 1 / (1 - var_total)

  # Return results as a list
  return(list(
    theta = theta,
    se_theta = se_theta,
    vcov_theta = vcov_theta,
    var_total = var_total,
    var_wtheta = var_wtheta,
    max_improvement_ratio = max_improvement_ratio
  ))
}