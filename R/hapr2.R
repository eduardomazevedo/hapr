#' HAPR second stage fit
#'
#' After fitting the first stage, we can specify an improvement ratio to estimate
#' the full model. Alternatively, we can specify the R-squared of the future fit
#' to estimate the improvement ratio. If specifying the r2_future, you can also
#' specify the r2_current to estimate the improvement ratio. Otherwise, the
#' r2_current is extracted from the data in the first stage fit.
#'
#' @param first_stage A hapr_first_stage_fit object
#' @param improvement_ratio The ratio to extrapolate by
#' @param r2_current The R-squared of the current fit
#' @param r2_future The R-squared of the future fit
#' @return A hapr_lm_fit object containing the results of the second stage
#' @export
hapr_second_stage <- function(
    first_stage,
    improvement_ratio = NULL,
    r2_current = NULL,
    r2_future = NULL,
    ...) {
  # Make sure first_stage is a hapr_lm_first_stage_fit object
  if (!inherits(first_stage, "hapr_first_stage_fit")) {
    stop("first_stage must be a hapr_first_stage_fit object.")
  }

  # Ensure that one and only one of improvement_ratio or r2_future is provided
  if (is.null(improvement_ratio) && is.null(r2_future)) {
    stop("Either improvement_ratio or r2_future must be specified.")
  }
  if (!is.null(improvement_ratio) && !is.null(r2_future)) {
    stop("Only one of improvement_ratio or r2_future should be provided.")
  }

  # If r2_current is not provided, extract it from first_stage
  if (is.null(r2_current)) {
    r2_current_source <- "first_stage"
    r2_current <- first_stage$regressions$y_on_gc$r2 |> as.numeric()
  } else {
    r2_current_source <- "user_provided"
  }

  # Compute the missing value
  if (is.null(r2_future)) {
    heritability_source <- "improvement_ratio"
    r2_future <- improvement_ratio * r2_current
  } else {
    heritability_source <- "r2_future"
    improvement_ratio <- r2_future / r2_current
  }

  if (improvement_ratio >= first_stage$stats$max_improvement_ratio) {
    stop(
      sprintf(
        "Improvement ratio must be less than %s.",
        first_stage$stats$max_improvement_ratio
      )
    )
  }

  # Var epsilon
  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- first_stage$stats$var_v_plus_var_epsilon - var_epsilon

  # Calculate a, b, and c
  posterior <- abc(var_epsilon, var_v)

  # beta calculation: collect coefficients
  gamma <- first_stage$coefficients$gamma
  theta <- first_stage$coefficients$theta
  beta <- gamma
  i_gc <- which(names(gamma) == "gc")
  i_other <- which(names(gamma) != "gc")

  # beta calculation: main logic
  if (first_stage$model_type == "lm") {
    beta[i_gc] <- gamma[i_gc] / posterior$a
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta
  } else if (first_stage$model_type == "probit") {
    sqrt_input <- posterior$a^2 - (gamma[i_gc]^2) * (posterior$c^2)
    if (sqrt_input < 0) {
        stop("Invalid posterior parameters: sqrt_input is negative")
        message(sprintf(
        "a = %f\ngamma_gc = %f\nc = %f\nsqrt_input = a^2 - (gamma_gc^2 * c^2) = %f",
        posterior$a,
        gamma[i_gc],
        posterior$c,
        sqrt_input
        ))
    }

    beta[i_gc] <-
        gamma[i_gc] / sqrt(sqrt_input)
    beta[i_other] <-
        gamma[i_other] * sqrt(1 + (posterior$c^2) * (beta[i_gc]^2)) -
        posterior$b * theta * beta[i_gc]
  }

  # Rename gc to gf in beta coefficients
  names(beta)[i_gc] <- "gf"

  # Update stripped model
  first_stage$regressions$y_on_gf_w$stripped_model$coefficients <- beta
  first_stage$regressions$y_on_gc_w$coefficients <- beta

  # Additional parameters
  if (first_stage$model_type == "lm") {
    additional_parameters <- list(
      var_eta = first_stage$regressions$y_on_gc_w$sigma_squared - posterior$c^2
    )
  } else if (first_stage$model_type == "probit") {
    additional_parameters <- list()
  }
  
  # Create the result object
  result <- list(
    model_type = first_stage$model_type,
    regressions = first_stage$regressions,
    coefficients = c(first_stage$coefficients, list(
      beta = beta
    )),
    additional_parameters = additional_parameters,
    stats = c(first_stage$stats, list(
      var_v = var_v,
      var_epsilon = var_epsilon,
      posterior = posterior,
      improvement_ratio = improvement_ratio,
      r2_current = r2_current,
      r2_future = r2_future,
      heritability_source = heritability_source,
      r2_current_source = r2_current_source
    ))
  )
  class(result) <- c("hapr_fit")
  result
}
