#' HAPR first stage fit
#'
#' @description
#' Fits the first stage of the HAPR model given the outcome y, PRS gc, and control variables w.
#' Returns Cov(gamma, theta) placeholder (currently NULL, can be extended with bootstrap).
#'
#' @param y Outcome variable
#' @param gc Polygenic risk score (normalized by preprocess)
#' @param w Data frame of control variables
#' @param model_type "lm", "probit", or "cox"
#' @return A hapr_first_stage_fit object
#' @export
hapr_first_stage <- function(y, gc, w, model_type) {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = model_type)
  y <- preprocessed$y; gc <- preprocessed$gc; w <- preprocessed$w

  # Define regression function by model type
  if (model_type == "lm") {
    regression_function <- function(data) strip_lm(stats::lm(y ~ ., data = data))
  } else if (model_type == "probit") {
    regression_function <- function(data) strip_probit(stats::glm(y ~ ., data = data,
                                                                  family = stats::binomial(link = "probit")))
  } else if (model_type == "cox") {
    regression_function <- function(data) {
      full_model <- survival::coxph(y ~ ., data = data)
      stripped <- strip_cox(full_model)
      stripped$model <- full_model  # keep full model for vcov()
      stripped
    }
  } else {
    stop("Unsupported model_type: ", model_type)
  }

  # Run regressions
  regressions <- list(
    gc_on_w   = strip_lm(stats::lm(gc ~ ., data = w)),
    y_on_w    = regression_function(w),
    y_on_gc   = regression_function(data.frame(gc = gc)),
    y_on_gc_w = regression_function(cbind(gc = gc, w)),
    y_on_gf_w = regression_function(cbind(gf = gc, w))
  )

  # Coefficients and vcovs
  coefficients <- list(
    theta = regressions$gc_on_w$coefficients,
    vcov_theta = regressions$gc_on_w$vcov_coefficients,
    gamma = regressions$y_on_gc_w$coefficients,
    vcov_gamma = {
      if (model_type == "cox") stats::vcov(regressions$y_on_gc_w$model)
      else regressions$y_on_gc_w$vcov_coefficients
    },
    vcov_gamma_theta = NULL  # placeholder for cross-covariance (can add later)
  )

  # Stats
  stats <- list(
    var_v_plus_var_epsilon = regressions$gc_on_w$sigma_squared,
    max_improvement_ratio  = 1 / (1 - regressions$gc_on_w$sigma_squared),
    var_wtheta             = regressions$gc_on_w$explained_variance
  )
  if (isTRUE(stats$var_v_plus_var_epsilon > 1)) {
    warning("The variance of v + epsilon is numerically > 1.")
    stats$var_v_plus_var_epsilon <- pmin(1 - stats$var_wtheta, 1)
    stats$max_improvement_ratio  <- Inf
  }

  result <- list(
    model_type = model_type,
    regressions = regressions,
    coefficients = coefficients,
    stats = stats
  )
  class(result) <- "hapr_first_stage_fit"
  result
}
