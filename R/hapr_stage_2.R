#' HAPR second stage fit
#'
#' @description
#' Fits the full HARP model given the first stage fit and an improvement ratio.
#' @importFrom numDeriv jacobian
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
    r2_future = NULL) {

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
    stop(sprintf("Improvement ratio must be less than %s.",
                 first_stage$stats$max_improvement_ratio))
  }

  derived_parameters <- calculate_parameters(first_stage$model_type, first_stage$coefficients, improvement_ratio)

  beta <- derived_parameters$beta
  var_v <- derived_parameters$var_v
  var_epsilon <- derived_parameters$var_epsilon
  posterior <- derived_parameters$posterior

  # --- Delta method for vcov and CI ---
  gamma_hat <- first_stage$coefficients$gamma
  theta_hat <- first_stage$coefficients$theta
  var_v_plus_var_epsilon_hat <- first_stage$coefficients$var_v_plus_var_epsilon
  vcov_gamma <- first_stage$vcov_coefficients$gamma
  vcov_theta <- first_stage$vcov_coefficients$theta
  vcov_var_v_plus_var_epsilon <- first_stage$vcov_coefficients$var_v_plus_var_epsilon

  # Store names
  names_gamma <- names(gamma_hat)
  names_theta <- names(theta_hat)
  stopifnot(all.equal(names_gamma, rownames(vcov_gamma)))
  stopifnot(all.equal(names_theta, rownames(vcov_theta)))
  stopifnot(all.equal(names_gamma, colnames(vcov_gamma)))
  stopifnot(all.equal(names_theta, colnames(vcov_theta)))

  param_hat <- c(gamma_hat, theta_hat, var_v_plus_var_epsilon_hat)
  names(param_hat) <- c(names(gamma_hat), names(theta_hat), "var_v_plus_var_epsilon")
  ng <- length(gamma_hat)
  nt <- length(theta_hat)

  vcov_full <- matrix(0, ng + nt + 1, ng + nt + 1)
  vcov_full[1:ng, 1:ng] <- vcov_gamma
  vcov_full[(ng + 1):(ng + nt), (ng + 1):(ng + nt)] <- vcov_theta
  vcov_full[(ng + nt + 1), (ng + nt + 1)] <- vcov_var_v_plus_var_epsilon
  rownames(vcov_full) <- colnames(vcov_full) <- names(param_hat)

  beta_wrapper <- function(params) {
    gamma <- params[1:ng]
    theta <- params[(ng + 1):(ng + nt)]
    var_v_plus_var_epsilon <- params[ng + nt + 1]
    names(gamma) <- names(gamma_hat)
    names(theta) <- names(theta_hat)
    derived_parameters <- calculate_parameters(first_stage$model_type, list(gamma = gamma, theta = theta, var_v_plus_var_epsilon = var_v_plus_var_epsilon), improvement_ratio)
    derived_parameters$beta
  }

  J_raw <- numDeriv::jacobian(beta_wrapper, param_hat)
  rownames(J_raw) <- names(beta_wrapper(param_hat))
  colnames(J_raw) <- names(param_hat)

  vcov_beta <- J_raw %*% vcov_full %*% t(J_raw)
  sd_beta <- sqrt(diag(vcov_beta))
  names(sd_beta) <- names(beta)

  # Confidence intervals
  z <- qnorm(0.975)
  ci_beta <- data.frame(
    Estimate = beta,
    Std.Error = sd_beta,
    Lower = beta - z * sd_beta,
    Upper = beta + z * sd_beta,
    row.names = names(beta),
    check.names = FALSE
  )

  # Create model for y_on_gf_w
  y_on_gf_w <- first_stage$regressions$y_on_gc_w

  # Update model coefficients
  y_on_gf_w$coefficients <- beta
  y_on_gf_w$vcov_coefficients <- vcov_beta
  y_on_gf_w$stripped_model$coefficients <- beta

  additional_parameters <- list()
  if (first_stage$model_type == "lm") {
    additional_parameters$var_eta <- first_stage$regressions$y_on_gc_w$sigma_squared - posterior$c^2
  } else if (first_stage$model_type == "cox") {
    theta_intercept <- first_stage$coefficients$theta[1]
    base_hazard_conversion_ratio <- exp(
      beta["gf"]^2 * posterior$c^2 / 2 + beta["gf"] * theta_intercept * posterior$b
    )

    baseline_hazard <- first_stage$regressions$y_on_gc_w$baseline_hazard
    baseline_hazard$hazard <- baseline_hazard$hazard / base_hazard_conversion_ratio

    y_on_gf_w$baseline_hazard <- baseline_hazard

    additional_parameters$base_hazard_conversion_ratio <- base_hazard_conversion_ratio
    additional_parameters$baseline_hazard <- baseline_hazard
  }

  # Clean up not needed parts of y_on_gf_w based on model type
  if (first_stage$model_type == "cox") {
    y_on_gf_w$stripped_model$var <- NULL
    y_on_gf_w$stripped_model$iter <- NULL
    y_on_gf_w$stripped_model$means <- NULL
    y_on_gf_w$stripped_model$method <- NULL
    y_on_gf_w$stripped_model$assign <- NULL
    y_on_gf_w$stripped_model$timefix <- NULL
    y_on_gf_w$stripped_model$formula <- NULL
    y_on_gf_w$stripped_model$xlevels <- NULL
    y_on_gf_w$stripped_model$contrasts <- NULL
    y_on_gf_w$stripped_model$formula <- NULL
  } else if (first_stage$model_type == "lm") {
    y_on_gf_w$stripped_model$rank <- NULL
    y_on_gf_w$stripped_model$assign <- NULL
    y_on_gf_w$stripped_model$qr <- NULL
    y_on_gf_w$stripped_model$df.residual <- NULL
    y_on_gf_w$stripped_model$contrasts <- NULL
    y_on_gf_w$stripped_model$xlevels <- NULL
    y_on_gf_w$stripped_model$call <- NULL
    y_on_gf_w$stripped_model$terms <- NULL
  } else if (first_stage$model_type == "probit") {
    y_on_gf_w$stripped_model$R <- NULL
    y_on_gf_w$stripped_model$rank <- NULL
    y_on_gf_w$stripped_model$deviance <- NULL
    y_on_gf_w$stripped_model$aic <- NULL
    y_on_gf_w$stripped_model$null.deviance <- NULL
    y_on_gf_w$stripped_model$iter <- NULL
    y_on_gf_w$stripped_model$df.residual <- NULL
    y_on_gf_w$stripped_model$df.null <- NULL
    y_on_gf_w$stripped_model$converged <- NULL
    y_on_gf_w$stripped_model$boundary <- NULL
    y_on_gf_w$stripped_model$call <- NULL
    y_on_gf_w$stripped_model$formula <- NULL
    y_on_gf_w$stripped_model$terms <- NULL
    y_on_gf_w$stripped_model$offset <- NULL
    y_on_gf_w$stripped_model$control <- NULL
    y_on_gf_w$stripped_model$contrasts <- NULL
    y_on_gf_w$stripped_model$xlevels <- NULL
    y_on_gf_w$stripped_model$qr <- NULL
  }

  result <- list(
    model_type = first_stage$model_type,
    regressions = c(
      first_stage$regressions,
      list(y_on_gf_w = y_on_gf_w)),
    coefficients = c(
      first_stage$coefficients,
      list(beta = beta)),
    vcov_coefficients = c(
      first_stage$vcov_coefficients,
      list(beta = vcov_beta)),
    standard_errors = sd_beta,
    ci_beta = ci_beta,
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
  class(result) <- "hapr_fit"
  result
}


#' Calculate beta coefficients based on model type
#'
#' @param model_type The type of model ("lm", "probit" or "cox")
#' @param coefficients A list containing gamma and theta coefficients
#' @param posterior The posterior values from the abc function
#' @return A named numeric vector of beta coefficients
calculate_parameters <- function(model_type, coefficients, improvement_ratio) {
  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- coefficients$var_v_plus_var_epsilon - var_epsilon
  posterior <- abc(var_epsilon, var_v)

  gamma <- coefficients$gamma
  theta <- coefficients$theta
  beta <- gamma  # default fallback

  i_gc <- which(names(gamma) == "gc")
  i_other <- which(names(gamma) != "gc")

  if (model_type == "lm") {
    beta[i_gc] <- gamma[i_gc] / posterior$a
    stopifnot(all(names(theta) == names(gamma[i_other])))  # ensure alignment
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta

  } else if (model_type == "probit") {
    sqrt_input <- posterior$a^2 - (gamma[i_gc]^2) * (posterior$c^2)
    if (sqrt_input < 0) stop("Invalid posterior parameters: sqrt_input is negative")
    beta[i_gc] <- gamma[i_gc] / sqrt(sqrt_input)
    beta[i_other] <- gamma[i_other] * sqrt(1 + (posterior$c^2) * beta[i_gc]^2) -
      posterior$b * theta * beta[i_gc]

  } else if (model_type == "cox") {
    theta_intercept <- theta[1]
    theta_others <- theta[-1]
    beta[i_gc] <- gamma[i_gc] / posterior$a
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta_others
  }

  names(beta)[i_gc] <- "gf"
  
  list(
    var_epsilon = var_epsilon,
    var_v = var_v,
    posterior = posterior,
    beta = beta)
}
