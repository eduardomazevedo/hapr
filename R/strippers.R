#' Strip a linear model to essential components
#'
#' @description
#' This internal function takes a linear model object and extracts the essential
#' components while removing memory-intensive parts of the model.
#'
#' @param fit An object of class 'lm' (linear model)
#'
#' @return A list containing:
#' \itemize{
#'   \item coefficients: Model coefficients
#'   \item vcov_coefficients: Variance-covariance matrix of coefficients
#'   \item r2: R-squared value of the model
#'   \item explained_variance: Variance of fitted values
#'   \item sigma_squared: Residual variance
#'   \item var_outcome: Variance of the outcome variable
#'   \item stripped_model: A trimmed version of the original model object
#' }
#'
#' @keywords internal
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
  fit$qr <- NULL
  fit$weights <- NULL
  fit$data <- NULL

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


#' Strip a probit model to essential components
#'
#' @description
#' This internal function takes a probit GLM object and extracts the essential
#' components while removing memory-intensive parts of the model.
#'
#' @param fit An object of class 'glm' with probit link function
#'
#' @return A list containing:
#' \itemize{
#'   \item coefficients: Model coefficients
#'   \item vcov_coefficients: Variance-covariance matrix of coefficients
#'   \item r2: R-squared value on the liability scale
#'   \item stripped_model: A trimmed version of the original model object
#' }
#'
#' @keywords internal
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
  fit$qr <- NULL
  fit$linear.predictors <- NULL
  fit$data <- NULL
  fit$prior.weights <- NULL

  # Return simplified object
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    r2 = r2,
    stripped_model = fit
  )
}


#' Strip a Cox proportional hazards model to essential components
#'
#' @description
#' This internal function takes a Cox proportional hazards model object and extracts
#' the essential components while removing memory-intensive parts of the model.
#'
#' @param fit An object of class 'coxph' from the survival package
#'
#' @return A list containing:
#' \itemize{
#'   \item coefficients: Model coefficients
#'   \item vcov_coefficients: Variance-covariance matrix of coefficients
#'   \item baseline_hazard: Baseline hazard function
#'   \item stripped_model: A trimmed version of the original model object
#' }
#'
#' @keywords internal
strip_cox <- function(fit) {
  # Check if model is a Cox proportional hazards model
  if (!inherits(fit, "coxph")) {
    stop("Model must be a coxph object from the survival package.")
  }

  # Extract coefficients
  coefficients <- fit$coefficients

  # Extract variance-covariance matrix of coefficients
  vcov_coefficients <- vcov(fit)

  # Compute baseline hazard function
  baseline_hazard <- survival::basehaz(fit, centered = FALSE)

  # Strip unnecessary elements
  fit$model <- NULL
  fit$y <- NULL
  fit$x <- NULL
  fit$residuals <- NULL
  fit$linear.predictors <- NULL
  fit$fitted.values <- NULL
  fit$effects <- NULL
  fit$weights <- NULL
  fit$call <- NULL
  fit$terms <- NULL
  # fit$formula <- NULL
  # fit$means <- NULL
  fit$concordance <- NULL
  fit$loglik <- NULL
  fit$wald.test <- NULL
  fit$score <- NULL
  fit$rscore <- NULL
  fit$n <- NULL
  fit$nevent <- NULL

  # Return everything in a clean list
  list(
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    baseline_hazard = baseline_hazard,
    stripped_model = fit
  )
}
