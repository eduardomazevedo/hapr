#' HAPR first stage fit
#'
#' @description
#' Fits the first stage of the HARP model given the outcome y, PRS gc, and control variables w.
#' 
#' @details
#' This returns a first stage fit, which does not need to assume an improvement ratio. Run
#' hapr_second_stage(first_stage_fit, improvement_ratio) to specify an improvement ratio and get the full model.
#'
#' @param y Outcome variable
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#' @param model_type "lm", "probit", or "cox"
#'
#' @return A hapr_first_stage_fit object containing the results of the first stage.
#' @export
hapr_first_stage <- function(y, gc, w, model_type) {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = model_type)
  y <- preprocessed$y
  gc <- preprocessed$gc
  w <- preprocessed$w
  rm(preprocessed)
  
  # Define regression function by model type
  if (model_type == "lm") {
    regression_function <- function(data) fit_lm_lowlevel(y, data)
  } else if (model_type == "probit") {
    regression_function <- function(data) fit_probit_lowlevel(y, data)
  } else if (model_type == "cox") {
    regression_function <- function(data) fit_cox_lowlevel(y, data)
  } else {
    stop("Unsupported model_type: ", model_type)
  }
  
  # Run regressions using low-level functions
  regressions <- list(
    gc_on_w = fit_lm_lowlevel(gc, w),
    y_on_w = regression_function(w),
    y_on_gc = regression_function(data.frame(gc = gc)),
    y_on_gc_w = regression_function(cbind(gc = gc, w))
  )
  
  # Extract coefficients
  coefficients <- list(
    theta = regressions$gc_on_w$coefficients,
    var_v_plus_var_epsilon = regressions$gc_on_w$sigma_squared,
    gamma = regressions$y_on_gc_w$coefficients
  )

  # Calculate vcov of var_v_plus_var_epsilon
  degrees_of_freedom <- regressions$gc_on_w$df_residual
  v_cov_var_v_plus_var_epsilon <- 2 * regressions$gc_on_w$sigma_squared^2 / degrees_of_freedom

  vcov_coefficients <- list(
    theta = regressions$gc_on_w$vcov_coefficients,
    var_v_plus_var_epsilon = v_cov_var_v_plus_var_epsilon,
    gamma = regressions$y_on_gc_w$vcov_coefficients
  )
  
  # Summary statistics
  stats <- list(
    max_improvement_ratio = 1 / (1 - regressions$gc_on_w$sigma_squared),
    var_wtheta = regressions$gc_on_w$explained_variance
  )
  if (coefficients$var_v_plus_var_epsilon > 1) {
    warning("The variance of v plus epsilon is numerically greater than 1.")
    coefficients$var_v_plus_var_epsilon <- pmin(1 - stats$var_wtheta, 1)
    stats$max_improvement_ratio <- Inf
  }
  
  # Return
  result <- list(
    model_type = model_type,
    regressions = regressions,
    coefficients = coefficients,
    vcov_coefficients = vcov_coefficients,
    stats = stats
  )
  class(result) <- "hapr_first_stage_fit"
  result
}