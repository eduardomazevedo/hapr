#' HAPR Cox proportional hazards model second stage fit
#'
#' After fitting the first stage, we can specify an improvement ratio to estimate 
#' the full model. Alternatively, we can specify the R-squared of the future fit 
#' to estimate the improvement ratio. If specifying the r2_future, you can also 
#' specify the r2_current to estimate the improvement ratio. Otherwise, the 
#' r2_current is extracted from the data in the first stage fit.
#'
#' @param first_stage A hapr_cox_first_stage_fit object
#' @param improvement_ratio The ratio to extrapolate by
#' @param r2_current The R-squared of the current fit
#' @param r2_future The R-squared of the future fit
#' @return A hapr_cox_fit object containing the results of the second stage
#' @export
hapr_cox_second_stage <- function(
  first_stage,
  improvement_ratio = NULL,
  r2_current = NULL,
  r2_future = NULL
) {
  # Make sure first_stage is a hapr_cox_first_stage_fit object
  if (!inherits(first_stage, "hapr_cox_first_stage_fit")) {
    stop("first_stage must be a hapr_cox_first_stage_fit object.")
  }

  # Ensure that one and only one of improvement_ratio or r2_future is provided
  if (is.null(improvement_ratio)) {
    if (is.null(r2_future) || is.null(r2_current)) {
      stop("Either improvement_ratio or r2_future and r2_current must be provided.")
    }
    improvement_ratio <- r2_future / r2_current
  }

  if (!is.null(improvement_ratio) && !is.null(r2_future)) {
    stop("Only one of improvement_ratio or r2_future should be provided.")
  }

  if (improvement_ratio >= first_stage$max_improvement_ratio) {
    stop(
      sprintf(
        "Improvement ratio must be less than %s.",
        first_stage$max_improvement_ratio
      )
    )
  }

  # Var epsilon
  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- first_stage$var_total - var_epsilon

  # Calculate a, b, and c
  posterior <- abc(var_epsilon, var_v)

  # beta
  gamma <- first_stage$gamma
  theta <- first_stage$theta
  theta_without_intercept <- theta[-1]
  beta <- gamma
  i_gc <- which(names(gamma) == 'gc')
  i_other <- which(names(gamma) != 'gc')
  beta[i_gc] <- gamma[i_gc] / posterior$a
  beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta_without_intercept

  # Rename gc to gf in beta coefficients
  names(beta)[i_gc] <- "gf"

  # Create the result object
  result <- c(first_stage, list(
    improvement_ratio = improvement_ratio,
    posterior_parameters = posterior,
    beta = beta
  ))

  class(result) <- "hapr_cox_fit"
  result
}