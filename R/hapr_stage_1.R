#' HAPR first stage fit
#'
#' @description
#' Fits the first stage of the HARP model given the outcome y, PRS gc, and control variables w.
#' 
#' @details
#' This returns a first stage fit, which does not need to assume an improvement ratio. Run
#' hapr_second_stage(first_stage_fit, improvement_ratio) to specify an improvement ratio and get the full model.
#' For model_type "mle", only the gc ~ w regression is run (used internally by
#' hapr_mle_survival).
#' 
#' **Coefficient Ordering:**
#' The coefficients are returned in a specific order that stage 2 relies on:
#' - `theta`: (Intercept), w1, w2, ... (from gc ~ w regression)
#' - `gamma`: gc, (Intercept), w1, w2, ... (from y ~ gc + w regression)
#' This ordering is critical for the stage 2 calculations.
#'
#' @param y Outcome variable. For "lm": numeric vector. For "probit": logical vector.
#' @param gc Polygenic risk score (numeric vector, will be normalized)
#' @param w Control variables (numeric matrix, must not include constant or linearly dependent columns)
#' @param model_type "lm", "probit", or "mle" (used internally by hapr_mle_survival)
#'
#' @return A hapr_first_stage_fit object containing the results of the first stage.
#' @export
hapr_first_stage <- function(y, gc, w, model_type) {
  # Validate model_type
  if (!model_type %in% c("lm", "probit", "mle")) {
    stop("model_type must be one of: 'lm', 'probit', 'mle'")
  }
  
  # Validate gc
  if (!is.numeric(gc)) {
    stop("gc must be numeric")
  }
  gc <- as.numeric(gc)
  
  # Validate w
  if (!is.matrix(w)) {
    stop("w must be a numeric matrix")
  }
  if (ncol(w) == 0) {
    stop("w must have at least one column")
  }
  if (!is.numeric(w)) {
    stop("w must be a numeric matrix")
  }
  
  if (any(is.na(gc))) {
    stop("gc contains missing values")
  }
  if (any(is.na(w))) {
    stop("w contains missing values")
  }
  
  # Check dimensions
  n <- length(gc)
  if (nrow(w) != n) {
    stop("gc and w must have the same number of observations")
  }
  
  # Validate y based on model_type
  if (model_type == "lm") {
    if (!is.numeric(y)) {
      stop("For 'lm' model_type, y must be numeric")
    }
    y <- as.numeric(y)
  } else if (model_type == "probit") {
    if (!is.logical(y)) {
      stop("For 'probit' model_type, y must be a logical vector")
    }
    # Convert logical to numeric for glm.fit (0/1)
    y <- as.numeric(y)
  } else {
    if (!missing(y) && length(y) != n) {
      stop("y, gc, and w must have the same number of observations")
    }
  }
  
  if (model_type != "mle") {
    if (length(y) != n) {
      stop("y, gc, and w must have the same number of observations")
    }
    if (any(is.na(y))) {
      stop("y contains missing values")
    }
  }
  
  # Check for constant columns in w
  w_var <- apply(w, 2, var, na.rm = TRUE)
  if (any(w_var == 0 | is.na(w_var))) {
    stop("w contains constant columns (zero variance). Remove constant columns before calling hapr_first_stage.")
  }
  
  # Check for linearly dependent columns in w (check rank)
  if (ncol(w) > 1) {
    w_rank <- qr(w)$rank
    if (w_rank < ncol(w)) {
      stop("w contains linearly dependent columns. Remove linearly dependent columns before calling hapr_first_stage.")
    }
  }
  
  # Always set column names to w1, w2, etc.
  colnames(w) <- paste0("w", seq_len(ncol(w)))
  
  # Normalize gc (scale to mean 0, variance 1)
  # gc <- as.numeric(scale(gc))
  
  # Helper function to create design matrix with named intercept
  add_intercept <- function(X) {
    X_with_int <- cbind(1, X)
    colnames(X_with_int)[1] <- "(Intercept)"
    X_with_int
  }
  
  # Run regressions using low-level functions
  # Create design matrices with named intercept column
  gc_on_w_X <- add_intercept(w)
  if (model_type == "mle") {
    regressions <- list(
      gc_on_w = fit_lm(gc, gc_on_w_X)
    )
  } else {
    # Define regression function by model type
    if (model_type == "lm") {
      regression_function <- function(X_mat) fit_lm(y, X_mat)
    } else if (model_type == "probit") {
      regression_function <- function(X_mat) fit_probit(y, X_mat)
    } else {
      stop("Unsupported model_type: ", model_type)
    }
    
    y_on_w_X <- add_intercept(w)
    y_on_gc_X <- add_intercept(gc)
    y_on_gc_w_X <- cbind(gc, add_intercept(w))
    colnames(y_on_gc_w_X)[1] <- "gc"
    
    regressions <- list(
      gc_on_w = fit_lm(gc, gc_on_w_X),
      y_on_w = regression_function(y_on_w_X),
      y_on_gc = regression_function(y_on_gc_X),
      y_on_gc_w = regression_function(y_on_gc_w_X)
    )
  }
  
  # Extract coefficients
  if (model_type == "mle") {
    parameters <- list(
      theta = regressions$gc_on_w$coefficients,
      var_v_plus_var_epsilon = regressions$gc_on_w$sigma_squared,
      gamma = NULL
    )
  } else {
    parameters <- list(
      theta = regressions$gc_on_w$coefficients,
      var_v_plus_var_epsilon = regressions$gc_on_w$sigma_squared,
      gamma = regressions$y_on_gc_w$coefficients
    )
  }

  # Extract vcov of var_v_plus_var_epsilon (now returned directly from fit_lm)
  v_cov_var_v_plus_var_epsilon <- regressions$gc_on_w$var_sigma_squared

  if (model_type == "mle") {
    vcov_parameters <- list(
      theta = regressions$gc_on_w$vcov_coefficients,
      var_v_plus_var_epsilon = v_cov_var_v_plus_var_epsilon,
      gamma = NULL
    )
  } else {
    vcov_parameters <- list(
      theta = regressions$gc_on_w$vcov_coefficients,
      var_v_plus_var_epsilon = v_cov_var_v_plus_var_epsilon,
      gamma = regressions$y_on_gc_w$vcov_coefficients
    )
  }
  
  # Calculate explained variance (variance of fitted values) for var_wtheta
  # For gc_on_w: explained_variance = r2 * var(gc)
  explained_variance <- regressions$gc_on_w$r2
  
  # Summary statistics
  stats <- list(
    max_improvement_ratio = 1 / (1 - regressions$gc_on_w$sigma_squared),
    var_wtheta = explained_variance
  )
  if (parameters$var_v_plus_var_epsilon > 1) {
    warning("The variance of v plus epsilon is numerically greater than 1.")
    parameters$var_v_plus_var_epsilon <- pmin(1 - stats$var_wtheta, 1)
    stats$max_improvement_ratio <- Inf
  }
  
  # Return
  result <- list(
    model_type = model_type,
    regressions = regressions,
    parameters = parameters,
    vcov_parameters = vcov_parameters,
    stats = stats,
    preprocessed = list(y = y, gc = gc, w = w)
  )
  class(result) <- "hapr_first_stage_fit"
  result
}
