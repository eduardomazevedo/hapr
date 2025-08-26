#' HAPR second stage fit
#'
#' @description
#' Fits the full HAPR model given the first stage fit and an improvement ratio or r2_future.
#'
#' @param first_stage A hapr_first_stage_fit object
#' @param improvement_ratio Ratio to extrapolate by (optional if r2_future supplied)
#' @param r2_current Optional R² of the current fit
#' @param r2_future Optional R² of the future fit (implies improvement_ratio)
#' @return A hapr_fit object
#' @export
hapr_second_stage <- function(
    first_stage,
    improvement_ratio = NULL,
    r2_current = NULL,
    r2_future = NULL
) {
  if (!inherits(first_stage, "hapr_first_stage_fit")) {
    stop("first_stage must be a hapr_first_stage_fit object.")
  }
  if (is.null(improvement_ratio) && is.null(r2_future)) {
    stop("Either improvement_ratio or r2_future must be specified.")
  }
  if (!is.null(improvement_ratio) && !is.null(r2_future)) {
    stop("Only one of improvement_ratio or r2_future should be provided.")
  }

  # R² bookkeeping
  if (is.null(r2_current)) {
    r2_current <- first_stage$regressions$y_on_gc$r2 |> as.numeric()
    if (first_stage$model_type == "cox") r2_current <- NA
  }
  if (is.null(r2_future)) {
    r2_future <- if (first_stage$model_type == "cox") NA else improvement_ratio * r2_current
  } else {
    improvement_ratio <- r2_future / r2_current
  }

  if (!is.na(improvement_ratio) &&
      improvement_ratio >= first_stage$stats$max_improvement_ratio) {
    stop("Improvement ratio too large")
  }

  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- first_stage$stats$var_v_plus_var_epsilon - var_epsilon
  posterior <- abc(var_epsilon, var_v)

  # Beta estimates
  beta <- calculate_beta(first_stage$model_type, first_stage$coefficients, posterior)

  # --- Delta method variance ---
  gamma_hat <- first_stage$coefficients$gamma
  theta_hat <- first_stage$coefficients$theta
  vcov_gamma <- first_stage$coefficients$vcov_gamma
  vcov_theta <- first_stage$coefficients$vcov_theta

  ng <- length(gamma_hat); nt <- length(theta_hat)
  vcov_full <- matrix(0, ng+nt, ng+nt)
  vcov_full[1:ng, 1:ng] <- vcov_gamma
  vcov_full[(ng+1):(ng+nt), (ng+1):(ng+nt)] <- vcov_theta

  param_hat <- c(gamma_hat, theta_hat)
  names(param_hat) <- c(names(gamma_hat), names(theta_hat))

  beta_wrapper <- function(params) {
    gamma <- params[1:ng]; theta <- params[(ng+1):(ng+nt)]
    names(gamma) <- names(gamma_hat); names(theta) <- names(theta_hat)
    calculate_beta(first_stage$model_type, list(gamma=gamma, theta=theta), posterior)
  }

  J <- numDeriv::jacobian(beta_wrapper, param_hat)
  vcov_beta <- J %*% vcov_full %*% t(J)
  sd_beta <- sqrt(diag(vcov_beta))
  names(sd_beta) <- names(beta)

  z <- stats::qnorm(0.975)
  ci_beta <- data.frame(
    Estimate = beta,
    Std.Error = sd_beta,
    Lower = beta - z*sd_beta,
    Upper = beta + z*sd_beta,
    row.names = names(beta),
    check.names = FALSE
  )

  result <- list(
    model_type = first_stage$model_type,
    regressions = first_stage$regressions,
    coefficients = c(first_stage$coefficients, list(beta = beta)),
    standard_errors = sd_beta,
    vcov_beta = vcov_beta,
    ci_beta = ci_beta,
    stats = c(first_stage$stats, list(
      var_v = var_v,
      var_epsilon = var_epsilon,
      posterior = posterior,
      improvement_ratio = improvement_ratio,
      r2_current = r2_current,
      r2_future = r2_future
    ))
  )
  class(result) <- "hapr_fit"
  result
}

#' Calculate beta coefficients based on model type
#'
#' Internal helper used by hapr_second_stage().
#' @keywords internal
calculate_beta <- function(model_type, coefficients, posterior) {
  gamma <- coefficients$gamma
  theta <- coefficients$theta
  beta <- gamma

  i_gc <- which(names(gamma) == "gc")
  i_other <- which(names(gamma) != "gc")

  if (model_type == "lm") {
    beta[i_gc] <- gamma[i_gc] / posterior$a
    stopifnot(all(names(theta) == names(gamma[i_other])))
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta

  } else if (model_type == "probit") {
    sqrt_input <- posterior$a^2 - (gamma[i_gc]^2) * (posterior$c^2)
    if (sqrt_input < 0) stop("Invalid posterior parameters: sqrt_input is negative")
    beta[i_gc] <- gamma[i_gc] / sqrt(sqrt_input)
    beta[i_other] <- gamma[i_other] * sqrt(1 + (posterior$c^2) * beta[i_gc]^2) -
      posterior$b * theta * beta[i_gc]

  } else if (model_type == "cox") {
    theta_others <- theta[-1]
    beta[i_gc] <- gamma[i_gc] / posterior$a
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta_others
  }

  names(beta)[i_gc] <- "gf"
  beta[order(names(beta))]
}
