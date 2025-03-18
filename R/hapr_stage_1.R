#' HAPR first stage fit
#' @param y Outcome variable
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#' @param model_type "lm" or "probit"
#'
#' @return A hapr_first_stage_fit object containing the results of the first stage.
#' @details
#' Fits the HARP model given the outcome y, PRS gc, and control variables w. This returns
#' a first stage fit, which does not need to assume an improvement ratio. Run
#' hapr_lm_second_stage(first_stage_fit, improvement_ratio) to specify an improvement ratio and get the full model.
#' @export
hapr_first_stage <- function(y, gc, w, model_type) {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = model_type)
  y <- preprocessed$y
  gc <- preprocessed$gc
  w <- preprocessed$w
  rm(preprocessed)

  if (model_type == "lm") {
    regression_function <- function(data) {
      strip_lm(lm(y ~ ., data = data))
    }
  } else if (model_type == "probit") {
    regression_function <- function(data) {
      strip_probit(glm(y ~ ., data = data, family = binomial(link = "probit")))
    }
  } else if (model_type == "cox") {
    regression_function <- function(data) {
      strip_cox(survival::coxph(y ~ ., data = data))
    }
  }

  # Regressions
  regressions <- list(
    gc_on_w = strip_lm(lm(gc ~ ., data = w)),
    y_on_w = regression_function(w),
    y_on_gc = regression_function(data.frame(gc = gc)),
    y_on_gc_w = regression_function(cbind(gc = gc, w)),
    y_on_gf_w = regression_function(cbind(gf = gc, w))
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
    model_type = model_type,
    regressions = regressions,
    coefficients = coefficients,
    stats = stats
  )
  class(result) <- "hapr_first_stage_fit"
  result
}
