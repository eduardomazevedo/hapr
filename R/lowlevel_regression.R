#' Low-level regression functions that return only essential data
#'
#' These functions use low-level R functions to fit models and return only
#' coefficients, variance-covariance matrices, terms objects (for prediction),
#' and necessary statistics, avoiding the bloat of full model objects.
#'
#' @noRd
NULL


#' Fit linear model on a matrix without saving bloat.
#'
#' @param y Numeric response vector
#' @param X Numeric matrix of predictors
#' @return List of model statistics including coefficients, vcov_coefficients, sigma_squared, var_sigma_squared, r2
#' @keywords internal
fit_lm <- function(y, X) {
  # 2. Fit model using the fastest base function
  fit <- .lm.fit(X, y)

  # 3. Dimensions
  n <- length(y)
  p <- ncol(X)
  df_residual <- n - p
  
  # 4. Calculate Sigma Squared (Residual Variance)
  residuals <- fit$residuals
  rss <- c(crossprod(residuals))
  sigma_squared <- rss / df_residual
  
  # 5. Calculate Variance of the Sigma Squared Estimator
  # Var(s^2) = 2 * sigma^4 / df
  var_sigma_squared <- 2 * (sigma_squared^2) / df_residual
  
  # 6. Coefficient VCOV: sigma^2 * (X'X)^-1
  # Handle potential LAPACK pivoting
  piv <- fit$pivot[1:p]
  
  XtX_inv_pivoted <- chol2inv(fit$qr[1:p, 1:p, drop = FALSE])
  
  XtX_inv <- matrix(0, p, p)
  XtX_inv[piv, piv] <- XtX_inv_pivoted
  
  vcov_coefficients <- XtX_inv * sigma_squared
  
  # 7. Final Formatting
  coef_names <- colnames(X)
  if (is.null(coef_names)) coef_names <- paste0("x", 1:p)
  
  dimnames(vcov_coefficients) <- list(coef_names, coef_names)
  coefficients <- fit$coefficients
  names(coefficients) <- coef_names
  
  # R-squared
  tss <- c(crossprod(y - mean(y)))
  r2 <- 1 - (rss / tss)

  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    sigma_squared = sigma_squared,
    var_sigma_squared = var_sigma_squared,
    r2 = r2
  )
}

#' Fit probit model on a matrix without saving bloat.
#'
#' @param y Binary response vector (numeric 0/1 or logical, will be converted to numeric)
#' @param X Numeric matrix of predictors (should include intercept column if needed)
#' @return List of model statistics including coefficients, vcov_coefficients, r2
#' @keywords internal
fit_probit <- function(y, X) {
  
  # 1. Fit model using glm.fit (faster than glm, bypasses formula parsing)
  # glm.fit takes x and y directly.
  fit <- glm.fit(x = X, y = y, family = binomial(link = "probit"))
  
  # Check for convergence
  if (!fit$converged) {
    warning("glm.fit: algorithm did not converge")
  }

  # 3. Extract Coefficients and Covariance
  # In GLM, vcov = (X'WX)^-1 * phi. For binomial, phi=1.
  # We calculate it directly from the QR decomposition of the weighted design matrix
  # which glm.fit stores in fit$qr.
  p <- ncol(X)
  piv <- fit$qr$pivot[1:p]
  
  XtWX_inv_pivoted <- chol2inv(fit$qr$qr[1:p, 1:p, drop = FALSE])
  
  XtWX_inv <- matrix(0, p, p)
  XtWX_inv[piv, piv] <- XtWX_inv_pivoted
  
  vcov_coefficients <- XtWX_inv # Dispersion parameter is 1 for binomial
  
  # 4. Final Formatting
  coef_names <- colnames(X)
  if (is.null(coef_names)) coef_names <- paste0("x", 1:p)
  
  dimnames(vcov_coefficients) <- list(coef_names, coef_names)
  coefficients <- fit$coefficients
  names(coefficients) <- coef_names
  
  # 5. Calculate Liability-Scale R-squared
  # Formula: Var(Xb) / (Var(Xb) + 1)
  # We isolate the linear predictor associated with predictors (excluding intercept)
  # Check if first column is intercept (by name or by being constant)
  has_intercept <- "(Intercept)" %in% coef_names || 
                   (p > 0 && length(unique(X[, 1])) == 1)
  
  if (has_intercept && p > 1) {
    # If intercept is present, remove it to calculate variance of explained part
    beta_predictors <- coefficients[-1]
    X_predictors <- X[, -1, drop = FALSE]
    
    if (length(beta_predictors) == 0) {
      r2 <- 0
    } else {
      # Calculate variance of the linear predictor (X * beta)
      # We can use var() directly on the vector
      linear_predictor <- c(X_predictors %*% beta_predictors)
      var_lp <- var(linear_predictor)
      r2 <- var_lp / (var_lp + 1)
    }
  } else {
    # If no intercept, all columns contribute to variance
    # Note: R2 definition without intercept is tricky, but this standardizes latent variance
    linear_predictor <- c(X %*% coefficients)
    var_lp <- var(linear_predictor)
    r2 <- var_lp / (var_lp + 1)
  }

  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    r2 = r2
  )
}
