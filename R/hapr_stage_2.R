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
  if (!inherits(first_stage, "hapr_first_stage_fit")) {
    stop("first_stage must be a hapr_first_stage_fit object.")
  }

  if (is.null(improvement_ratio) && is.null(r2_future)) {
    stop("Either improvement_ratio or r2_future must be specified.")
  }
  if (!is.null(improvement_ratio) && !is.null(r2_future)) {
    stop("Only one of improvement_ratio or r2_future should be provided.")
  }

  if (is.null(r2_current)) {
    r2_current_source <- "first_stage"
    r2_current <- first_stage$regressions$y_on_gc$r2 |> as.numeric()
    if (first_stage$model_type == "cox") {
      r2_current <- NA
      r2_current_source <- "cox model, not available"
    }
  } else {
    r2_current_source <- "user_provided"
  }

  if (is.null(r2_future)) {
    heritability_source <- "improvement_ratio"
    r2_future <- improvement_ratio * r2_current
    if (first_stage$model_type == "cox") {
      r2_future <- NA
    }
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

  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- first_stage$stats$var_v_plus_var_epsilon - var_epsilon
  posterior <- abc(var_epsilon, var_v)

  beta <- calculate_beta(first_stage$model_type, first_stage$coefficients, posterior)

  # Update stripped model
  first_stage$regressions$y_on_gf_w$stripped_model$coefficients <- beta
  first_stage$regressions$y_on_gc_w$coefficients <- beta

  # Compute additional parameters based on the model type
  additional_parameters <- list()

  if (first_stage$model_type == "lm") {
    additional_parameters$var_eta <- first_stage$regressions$y_on_gc_w$sigma_squared - posterior$c^2
  } else if (first_stage$model_type == "cox") {
    # Compute the base hazard conversion ratio
    theta_intercept <- first_stage$coefficients$theta[1]
    base_hazard_conversion_ratio <- exp(
      beta["gf"]^2 * posterior$c^2 / 2 +
        beta["gf"] * theta_intercept * posterior$b
    )

    # Adjust the baseline hazard
    baseline_hazard <- first_stage$regressions$y_on_gf_w$baseline_hazard
    baseline_hazard$hazard <- baseline_hazard$hazard / base_hazard_conversion_ratio

    # Update `first_stage` with the modified baseline hazard
    first_stage$regressions$y_on_gf_w$baseline_hazard <- baseline_hazard
    first_stage$regressions$y_on_gf_w$stripped_model$baseline_hazard <- baseline_hazard

    # Store results in `additional_parameters`
    additional_parameters$base_hazard_conversion_ratio <- base_hazard_conversion_ratio
    additional_parameters$baseline_hazard <- baseline_hazard
  }


  result <- list(
    model_type = first_stage$model_type,
    regressions = first_stage$regressions,
    coefficients = c(first_stage$coefficients, list(beta = beta)),
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

#' Calculate beta coefficients based on model type
#'
#' @param model_type The type of model ("lm" or "probit")
#' @param coefficients A list containing gamma and theta coefficients
#' @param posterior The posterior values from the abc function
#' @return A named numeric vector of beta coefficients
calculate_beta <- function(model_type, coefficients, posterior) {
  gamma <- coefficients$gamma
  theta <- coefficients$theta
  beta <- gamma

  i_gc <- which(names(gamma) == "gc")
  i_other <- which(names(gamma) != "gc")

  if (model_type == "lm") {
    beta[i_gc] <- gamma[i_gc] / posterior$a
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta
  } else if (model_type == "probit") {
    sqrt_input <- posterior$a^2 - (gamma[i_gc]^2) * (posterior$c^2)
    if (sqrt_input < 0) {
      stop("Invalid posterior parameters: sqrt_input is negative")
    }

    beta[i_gc] <- gamma[i_gc] / sqrt(sqrt_input)
    beta[i_other] <- gamma[i_other] * sqrt(1 + (posterior$c^2) * (beta[i_gc]^2)) -
      posterior$b * theta * beta[i_gc]
  } else if (model_type == "cox") {
    theta_intercept <- theta[1] # Remeber this will matter for the base hazard
    theta_without_intercept <- theta[-1]
    beta[i_gc] <- gamma[i_gc] / posterior$a
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta_without_intercept
  }

  # Rename gc to gf
  names(beta)[i_gc] <- "gf"

  beta
}
