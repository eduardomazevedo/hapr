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

#' Fit linear model for scaled gc with variance-constrained uncertainty.
#'
#' Uses OLS coefficients for gc ~ w, but enforces Var(gc)=1 by setting
#' sigma^2 = 1 - Var(w * theta_w), and propagates uncertainty jointly for
#' (theta, sigma^2) via the delta method.
#'
#' @param y Numeric response vector (scaled gc)
#' @param X Numeric design matrix including intercept
#' @param slope_idx Integer indices of slope columns in X (exclude intercept)
#' @return Same structure as fit_lm, plus joint covariance pieces for
#'   coefficients and sigma_squared.
#' @keywords internal
fit_lm_scaled_gc <- function(y, X, slope_idx) {
  base_fit <- fit_lm(y = y, X = X)

  p <- ncol(X)
  coef_names <- names(base_fit$coefficients)

  if (length(slope_idx) == 0) {
    stop("slope_idx must contain at least one non-intercept column index.")
  }
  if (any(slope_idx < 1 | slope_idx > p)) {
    stop("slope_idx contains invalid column indices.")
  }

  # Recover (X'X)^-1 from the OLS vcov and unscaled sigma^2.
  xtx_inv <- base_fit$vcov_coefficients / base_fit$sigma_squared

  X_slope <- X[, slope_idx, drop = FALSE]
  beta_slope <- base_fit$coefficients[slope_idx]
  sigma_w <- stats::cov(X_slope)
  explained_var <- as.numeric(t(beta_slope) %*% sigma_w %*% beta_slope)
  sigma_squared <- 1 - explained_var

  # Numerical guardrails for near-boundary designs.
  sigma_squared <- min(max(sigma_squared, .Machine$double.eps), 1 - .Machine$double.eps)

  vcov_coefficients <- xtx_inv * sigma_squared

  grad_sigma <- rep(0, p)
  grad_sigma[slope_idx] <- -2 * as.numeric(sigma_w %*% beta_slope)
  cov_coefficients_sigma_squared <- as.numeric(vcov_coefficients %*% grad_sigma)
  names(cov_coefficients_sigma_squared) <- coef_names
  var_sigma_squared <- as.numeric(crossprod(grad_sigma, vcov_coefficients %*% grad_sigma))

  vcov_joint <- matrix(0, nrow = p + 1, ncol = p + 1)
  vcov_joint[seq_len(p), seq_len(p)] <- vcov_coefficients
  vcov_joint[seq_len(p), p + 1] <- cov_coefficients_sigma_squared
  vcov_joint[p + 1, seq_len(p)] <- cov_coefficients_sigma_squared
  vcov_joint[p + 1, p + 1] <- var_sigma_squared
  dimnames(vcov_joint) <- list(c(coef_names, "sigma_squared"), c(coef_names, "sigma_squared"))

  base_fit$vcov_coefficients <- vcov_coefficients
  base_fit$sigma_squared <- sigma_squared
  base_fit$var_sigma_squared <- var_sigma_squared
  base_fit$r2 <- 1 - sigma_squared
  base_fit$cov_coefficients_sigma_squared <- cov_coefficients_sigma_squared
  base_fit$vcov_coefficients_sigma_squared <- vcov_joint

  base_fit
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
  # This is invariant to the intercept because adding a constant does not change variance.
  linear_predictor <- c(X %*% coefficients)
  var_lp <- var(linear_predictor)
  r2 <- var_lp / (var_lp + 1)

  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    r2 = r2
  )
}
