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
#' @param softmax_correction, only used if model_type is "cox". Can be "clt" (default), "softmax-fast", or "softmax-slow".
#'
#' @return A hapr_first_stage_fit object containing the results of the first stage.
#' @export
hapr_first_stage <- function(y, gc, w, model_type, softmax_correction = "clt") {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = model_type)
  y <- preprocessed$y
  gc <- preprocessed$gc
  w <- preprocessed$w
  rm(preprocessed)
  
  # Define regression function by model type
  if (model_type == "lm") {
    regression_function <- function(data) strip_lm(lm(y ~ ., data = data))
  } else if (model_type == "probit") {
    regression_function <- function(data) strip_probit(glm(y ~ ., data = data, family = binomial(link = "probit")))
  } else if (model_type == "cox") {
    regression_function <- function(data) strip_cox(survival::coxph(y ~ ., data = data), softmax_correction)
  } else {
    stop("Unsupported model_type: ", model_type)
  }

  # Validate softmax_correction
  if (model_type == "cox") {
    if (!(softmax_correction %in% c("clt", "softmax-fast", "softmax-slow"))) {
      stop("softmax_correction must be 'clt', 'softmax-fast', or 'softmax-slow' if model_type is 'cox'.")
    }
  }
  
  # Run regressions
  regressions <- list(
    gc_on_w = strip_lm(lm(gc ~ ., data = w)),
    y_on_w = regression_function(w),
    y_on_gc = regression_function(data.frame(gc = gc)),
    y_on_gc_w = regression_function(cbind(gc = gc, w))
  )
  
  # Extract coefficients
  coefficients <- list(
    theta = regressions$gc_on_w$coefficients,
    vcov_theta = regressions$gc_on_w$vcov_coefficients,
    gamma = regressions$y_on_gc_w$coefficients,
    vcov_gamma = {
      if (model_type == "cox") {
        # Use full Cox model for vcov
        regressions$y_on_gc_w$vcov_coefficients
      } else {
        regressions$y_on_gc_w$vcov_coefficients
      }
    }
  )

  if (model_type == "cox") {
    coefficients$psi_hat <- regressions$y_on_gc_w$psi_hat
  }
  
  # Summary statistics
  stats <- list(
    var_v_plus_var_epsilon = regressions$gc_on_w$sigma_squared,
    max_improvement_ratio = 1 / (1 - regressions$gc_on_w$sigma_squared),
    var_wtheta = regressions$gc_on_w$explained_variance
  )
  if (stats$var_v_plus_var_epsilon > 1) {
    warning("The variance of v plus epsilon is numerically greater than 1.")
    stats$var_v_plus_var_epsilon <- pmin(1 - stats$var_wtheta, 1)
    stats$max_improvement_ratio <- Inf
  }
  # Return
  result <- list(
    model_type = model_type,
    regressions = regressions,
    coefficients = coefficients,
    stats = stats
  )
  class(result) <- "hapr_first_stage_fit"
  result
}