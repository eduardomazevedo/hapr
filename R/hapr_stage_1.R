#' HAPR first stage fit
#' Fits the first stage of the HARP model given y, PRS gc, and controls w.
#' @export
hapr_first_stage <- function(y, gc, w, model_type) {
  preprocessed <- preprocess(y, gc, w, model_type = model_type)
  y  <- preprocessed$y
  gc <- preprocessed$gc
  w  <- preprocessed$w

  if (model_type == "lm") {
    regression_function <- function(data) strip_lm(stats::lm(y ~ ., data = data))
  } else if (model_type == "probit") {
    regression_function <- function(data) strip_probit(stats::glm(
      y ~ ., data = data, family = stats::binomial(link = "probit")
    ))
  } else if (model_type == "cox") {
    regression_function <- function(data) {
      full_model <- survival::coxph(y ~ ., data = data)
      stripped <- strip_cox(full_model)
      stripped$model <- full_model
      stripped
    }
  } else {
    stop("Unsupported model_type: ", model_type)
  }

  regressions <- list(
    gc_on_w    = strip_lm(stats::lm(gc ~ ., data = w)),
    y_on_w     = regression_function(w),
    y_on_gc    = regression_function(data.frame(gc = gc)),
    y_on_gc_w  = regression_function(cbind(gc = gc, w)),
    y_on_gf_w  = regression_function(cbind(gf = gc, w))
  )

  # extract and align
  gamma <- regressions$y_on_gc_w$coefficients
  theta <- regressions$gc_on_w$coefficients
  vcov_gamma <- if (model_type == "cox") stats::vcov(regressions$y_on_gc_w$model)
  else regressions$y_on_gc_w$vcov_coefficients
  vcov_theta <- regressions$gc_on_w$vcov_coefficients

  vcov_gamma <- vcov_gamma[match(names(gamma), rownames(vcov_gamma)),
                           match(names(gamma), colnames(vcov_gamma)), drop = FALSE]
  rownames(vcov_gamma) <- colnames(vcov_gamma) <- names(gamma)

  vcov_theta <- vcov_theta[match(names(theta), rownames(vcov_theta)),
                           match(names(theta), colnames(vcov_theta)), drop = FALSE]
  rownames(vcov_theta) <- colnames(vcov_theta) <- names(theta)

  coefficients <- list(
    gamma = gamma,
    theta = theta,
    vcov_gamma = vcov_gamma,
    vcov_theta = vcov_theta,
    vcov_gamma_theta = NULL
  )

  # Summary statistics (numerically stable, no warnings)
  v_raw <- regressions$gc_on_w$sigma_squared
  eps_cap  <- 1e-9
  v_clamped <- min(1 - eps_cap, max(0, v_raw))

  stats <- list(
    var_v_plus_var_epsilon = v_clamped,
    max_improvement_ratio  = 1 / (1 - v_clamped),
    var_wtheta             = regressions$gc_on_w$explained_variance
  )

  result <- list(
    model_type = model_type,
    regressions = regressions,
    coefficients = coefficients,
    stats = stats
  )
  class(result) <- "hapr_first_stage_fit"
  result
}
