#' HAPR linear model first stage fit
#' @param y Outcome variable
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#'
#' @return A hapr_lm_first_stage_fit object containing the results of the first stage.
#' @details
#' Fits the HARP model given the outcome y, PRS gc, and control variables w. This returns
#' a first stage fit, which does not need to assume an improvement ratio. Run
#' hapr_lm_second_stage(first_stage_fit, improvement_ratio) to specify an improvement ratio and get the full model.
#' @export
hapr_lm_first_stage <- function(y, gc, w) {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = "lm")
  y <- preprocessed$y
  gc <- preprocessed$gc
  w <- preprocessed$w
  rm(preprocessed)

  # Regressions
  regressions <- list(
    gc_on_w = strip_lm(lm(gc ~ ., data = w)),
    y_on_w = strip_lm(lm(y ~ ., data = w)),
    y_on_gc = strip_lm(lm(y ~ gc)),
    y_on_gc_w = strip_lm(lm(y ~ ., data = cbind(gc = gc, w)))
  )

  # First stage results
  coefficients <- list(
    theta = regressions$gc_on_w$coefficients,
    vcov_theta = regressions$gc_on_w$vcov_coefficients,
    gamma = regressions$y_on_gc_w$coefficients,
    vcov_gamma = regressions$y_on_gc_w$vcov_coefficients
  )
  stats <- list(
    var_v_plus_var_epsilon = regressions$gc_on_w$sigma_squared,
    max_improvement_ratio = 1 / (1 - regressions$gc_on_w$sigma_squared),
    var_wtheta = regressions$gc_on_w$explained_variance
  )

  # Return
  result <- list(
    regressions = regressions,
    coefficients = coefficients,
    stats = stats
  )
  class(result) <- "hapr_lm_first_stage_fit"
  result
}
