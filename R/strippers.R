strip_lm <- function(fit) {
  # Check if model is a linear model
  if (!inherits(fit, "lm")) {
    stop("Model must be a linear model.")
  }

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

  fit$model <- NULL
  fit$y <- NULL
  fit$x <- NULL
  fit$residuals <- NULL
  fit$fitted.values <- NULL
  fit$effects <- NULL
  fit$qr$qr <- NULL
  fit$weights <- NULL
  
  # Return everything in a clean list
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,    
    r2 = r2,
    explained_variance = explained_variance,
    sigma_squared = sigma_squared,
    var_outcome = outcome_variance,
    stripped_model = fit
  )
}


strip_probit <- function(fit) {
  # Check if model is a probit glm
  if (!inherits(fit, "glm") || fit$family$link != "probit") {
    stop("Model must be a glm with probit link.")
  }

  # Coefficients
  coefficients <- coef(fit)

  # Variance-covariance matrix of coefficients
  vcov_coefficients <- vcov(fit)

  # Liability-scale R²
  r2 <- r2_liability_probit(fit)

  # Strip unnecessary elements
  fit$model <- NULL
  fit$y <- NULL
  fit$x <- NULL
  fit$residuals <- NULL
  fit$fitted.values <- NULL
  fit$effects <- NULL
  fit$qr$qr <- NULL
  fit$weights <- NULL
  
  # Return simplified object
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    r2 = r2,
    stripped_model = fit
  )
}
