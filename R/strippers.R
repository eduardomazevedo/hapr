strip_lm <- function(fit) {
  # Coefficients
  coefficients <- fit$coefficients

  # Variance-covariance matrix of coefficients
  vcov_coefficients <- vcov(fit)

  # Compute R^2
  r2 <- summary(fit)$r.squared
  
  # Explained variance: variance of fitted values
  fitted_values <- fitted(fit)
  explained_variance <- var(fitted_values)
  
  # Residual variance (sigma^2)
  sigma_squared <- summary(fit)$sigma^2
  
  # Variance of outcome variable
  outcome_variance <- var(fitted_values + residuals(fit))
  
  # Stripped-down lm object
  strip_lm <- function(fit) {
    fit$model <- NULL
    fit$y <- NULL
    fit$x <- NULL
    fit$residuals <- NULL
    fit$fitted.values <- NULL
    fit$effects <- NULL
    fit$qr$qr <- NULL
    fit$weights <- NULL
    fit
  }
  
  fit_stripped <- strip_lm(fit)
  
  # Return everything in a clean list
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,    
    r2 = r2,
    explained_variance = explained_variance,
    sigma_squared = sigma_squared,
    var_outcome = outcome_variance,
    stripped_model = fit_stripped
  )
}
