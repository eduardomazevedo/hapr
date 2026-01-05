#' Low-level regression functions that return only essential data
#'
#' These functions use low-level R functions to fit models and return only
#' coefficients, variance-covariance matrices, terms objects (for prediction),
#' and necessary statistics, avoiding the bloat of full model objects.
#'
#' @noRd
NULL

#' Fit a linear model using low-level functions
#'
#' @param y Response vector
#' @param data Data frame of predictors
#' @return List with coefficients, vcov_coefficients, terms, sigma_squared, 
#'   explained_variance, df_residual, r2
#' @noRd
fit_lm_lowlevel <- function(y, data) {
  # Create model frame and terms for prediction
  mf <- model.frame(y ~ ., data = data)
  mt <- terms(mf)
  
  # Create design matrix (includes intercept)
  X <- model.matrix(mt, mf)
  
  # Fit using .lm.fit (fastest, most stable low-level function)
  fit <- .lm.fit(X, y)
  
  # Calculate degrees of freedom
  n <- length(y)
  p <- ncol(X)
  df_residual <- n - p
  
  # Calculate residuals
  residuals <- fit$residuals
  fitted_values <- y - residuals
  
  # Calculate sigma squared (residual variance)
  sigma_squared <- sum(residuals^2) / df_residual
  
  # Calculate explained variance (variance of fitted values)
  explained_variance <- var(fitted_values)
  
  # Calculate R-squared
  var_y <- var(y)
  r2 <- if (var_y > 0) explained_variance / var_y else 0
  
  # Calculate variance-covariance matrix of coefficients
  # vcov = sigma^2 * (X'X)^{-1}
  # Use QR decomposition for numerical stability
  qr_X <- qr(X)
  if (qr_X$rank < p) {
    # Rank-deficient case
    coef_names <- colnames(X)
    vcov_coefficients <- matrix(NA, p, p, dimnames = list(coef_names, coef_names))
    R_inv <- tryCatch(chol2inv(qr.R(qr_X)), error = function(e) solve(qr.R(qr_X)))
    vcov_coefficients[qr_X$pivot[1:qr_X$rank], qr_X$pivot[1:qr_X$rank]] <- 
      sigma_squared * R_inv
  } else {
    # Full rank case
    R_inv <- tryCatch(chol2inv(qr.R(qr_X)), error = function(e) solve(qr.R(qr_X)))
    vcov_coefficients <- sigma_squared * R_inv
    colnames(vcov_coefficients) <- rownames(vcov_coefficients) <- colnames(X)
  }
  
  # Get coefficients (may have NA for rank-deficient cases)
  coefficients <- fit$coefficients
  names(coefficients) <- colnames(X)
  
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    terms = mt,
    xlevels = .getXlevels(mt, mf),
    sigma_squared = sigma_squared,
    explained_variance = explained_variance,
    df_residual = df_residual,
    r2 = r2
  )
}

#' Fit a probit model (still needs glm for IRLS, but extracts only essentials)
#'
#' @param y Response vector (binary)
#' @param data Data frame of predictors
#' @return List with coefficients, vcov_coefficients, terms, r2
#' @noRd
fit_probit_lowlevel <- function(y, data) {
  # Create model frame and terms for prediction
  mf <- model.frame(y ~ ., data = data)
  mt <- terms(mf)
  
  # Fit using glm (no low-level alternative for IRLS)
  fit <- glm(y ~ ., data = data, family = binomial(link = "probit"))
  
  # Extract coefficients
  coefficients <- coef(fit)
  
  # Extract variance-covariance matrix
  vcov_coefficients <- vcov(fit)
  
  # Calculate liability-scale R²
  # Need to compute this from the design matrix and coefficients
  X <- model.matrix(mt, mf)[, -1, drop = FALSE] # Remove intercept
  beta_hat <- coefficients[-1] # Remove intercept
  
  if (length(beta_hat) == 0) {
    r2 <- 0
  } else {
    var_linear_predictor <- var(X %*% beta_hat)
    r2 <- var_linear_predictor / (var_linear_predictor + 1)
  }
  
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    terms = mt,
    xlevels = .getXlevels(mt, mf),
    r2 = r2
  )
}

#' Fit a Cox model (still needs coxph for basehaz, but extracts only essentials)
#'
#' @param y Surv object (response)
#' @param data Data frame of predictors
#' @return List with coefficients, vcov_coefficients, terms, baseline_hazard
#' @noRd
fit_cox_lowlevel <- function(y, data) {
  # Fit using coxph (required for basehaz)
  fit <- survival::coxph(y ~ ., data = data)
  
  # Extract coefficients
  coefficients <- coef(fit)
  
  # Extract variance-covariance matrix
  vcov_coefficients <- vcov(fit)
  
  # Extract baseline hazard (required for Cox models)
  baseline_hazard <- survival::basehaz(fit, centered = FALSE)
  
  # Extract terms object for prediction
  # coxph stores terms in the fit object
  mt <- fit$terms
  if (is.null(mt)) {
    # Fallback: create terms from the formula
    formula_obj <- formula(fit)
    mf <- model.frame(formula_obj, data = cbind(data.frame(.y = y), data))
    mt <- delete.response(terms(mf))
  }
  
  # Get xlevels from the fit object if available
  xlevels <- fit$xlevels
  if (is.null(xlevels)) {
    # Try to extract from model frame
    tryCatch({
      mf <- model.frame(fit)
      xlevels <- .getXlevels(mt, mf)
    }, error = function(e) {
      xlevels <- NULL
    })
  }
  
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    terms = mt,
    xlevels = xlevels,
    baseline_hazard = baseline_hazard
  )
}