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

  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- first_stage$stats$var_v_plus_var_epsilon - var_epsilon
  posterior <- abc(var_epsilon, var_v)

  # Handle alpha for Cox model
  alpha_hat <- NULL
  if (first_stage$model_type == "cox") {
    psi_hat <- first_stage$coefficients$psi_hat
    alpha_hat <- 1 - (posterior$c) ^ 2 * psi_hat / 2
  }

  beta <- calculate_beta(
    model_type = first_stage$model_type,
    coefficients = first_stage$coefficients,
    posterior = posterior,
    alpha = alpha_hat
  )

  # --- Delta method for vcov and CI ---
  gamma_hat <- first_stage$coefficients$gamma
  theta_hat <- first_stage$coefficients$theta
  vcov_gamma <- first_stage$coefficients$vcov_gamma
  vcov_theta <- first_stage$coefficients$vcov_theta

  gamma_hat <- gamma_hat[order(names(gamma_hat))]
  theta_hat <- theta_hat[order(names(theta_hat))]
  vcov_gamma <- vcov_gamma[order(rownames(vcov_gamma)), order(colnames(vcov_gamma))]
  vcov_theta <- vcov_theta[order(rownames(vcov_theta)), order(colnames(vcov_theta))]

  param_hat <- c(gamma_hat, theta_hat)
  names(param_hat) <- c(names(gamma_hat), names(theta_hat))
  ng <- length(gamma_hat)
  nt <- length(theta_hat)

  vcov_full <- matrix(0, ng + nt, ng + nt)
  vcov_full[1:ng, 1:ng] <- vcov_gamma
  vcov_full[(ng + 1):(ng + nt), (ng + 1):(ng + nt)] <- vcov_theta
  rownames(vcov_full) <- colnames(vcov_full) <- names(param_hat)

  beta_wrapper <- function(params) {
    gamma <- params[1:ng]
    theta <- params[(ng + 1):(ng + nt)]
    names(gamma) <- names(gamma_hat)
    names(theta) <- names(theta_hat)
    out <- calculate_beta(
      model_type = first_stage$model_type,
      coefficients = list(gamma = gamma, theta = theta, psi_hat = first_stage$coefficients$psi_hat),
      posterior = posterior,
      alpha = alpha_hat
    )
    out[order(names(out))]
  }

  J_raw <- numDeriv::jacobian(beta_wrapper, param_hat)
  rownames(J_raw) <- names(beta_wrapper(param_hat))
  colnames(J_raw) <- names(param_hat)

  J <- J_raw[match(names(beta), rownames(J_raw)),
             match(names(param_hat), colnames(J_raw))]
  vcov_ordered <- vcov_full[match(names(param_hat), rownames(vcov_full)),
                            match(names(param_hat), colnames(vcov_full))]

  vcov_beta <- J %*% vcov_ordered %*% t(J)
  sd_beta <- sqrt(diag(vcov_beta))
  names(sd_beta) <- names(beta)

  z <- qnorm(0.975)
  ci_beta <- data.frame(
    Estimate = beta,
    Std.Error = sd_beta,
    Lower = beta - z * sd_beta,
    Upper = beta + z * sd_beta,
    row.names = names(beta),
    check.names = FALSE
  )

  y_on_gf_w <- first_stage$regressions$y_on_gc_w
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

  if (first_stage$model_type == "cox") {
    y_on_gf_w$stripped_model[c("var", "iter", "means", "method", "assign",
                               "timefix", "formula", "xlevels", "contrasts")] <- NULL
  } else if (first_stage$model_type == "lm") {
    y_on_gf_w$stripped_model[c("rank", "assign", "qr", "df.residual",
                               "contrasts", "xlevels", "call", "terms")] <- NULL
  } else if (first_stage$model_type == "probit") {
    y_on_gf_w$stripped_model[c("R", "rank", "deviance", "aic", "null.deviance", "iter",
                               "df.residual", "df.null", "converged", "boundary", "call",
                               "formula", "terms", "offset", "control", "contrasts", "xlevels", "qr")] <- NULL
  }

  stats_list <- list(
    var_v = var_v,
    var_epsilon = var_epsilon,
    posterior = posterior,
    improvement_ratio = improvement_ratio,
    r2_current = r2_current,
    r2_future = r2_future,
    heritability_source = heritability_source,
    r2_current_source = r2_current_source
  )

  if (first_stage$model_type == "cox") {
    stats_list$psi_hat <- psi_hat
    stats_list$alpha_hat <- alpha_hat
  }

  result <- list(
    model_type = first_stage$model_type,
    regressions = c(first_stage$regressions, list(y_on_gf_w = y_on_gf_w)),
    coefficients = c(first_stage$coefficients, list(beta = beta)),
    standard_errors = sd_beta,
    vcov_beta = vcov_beta,
    ci_beta = ci_beta,
    additional_parameters = additional_parameters,
    stats = c(first_stage$stats, stats_list)
  )
  class(result) <- "hapr_fit"
  result
}

#' Calculate beta coefficients based on model type
#'
#' @param model_type The type of model ("lm", "probit" or "cox")
#' @param coefficients A list containing gamma and theta coefficients (and psi_hat for cox)
#' @param posterior The posterior values from the abc function
#' @param alpha Optional alpha softmax correction, is only passed for cox model.
#' @return A named numeric vector of beta coefficients
calculate_beta <- function(model_type, coefficients, posterior, alpha = NULL) {
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
    if (is.null(alpha)) stop("Alpha must be provided for Cox model.")
    gamma <- gamma / alpha

    theta_intercept <- theta[1]
    theta_others <- theta[-1]
    beta[i_gc] <- gamma[i_gc] / posterior$a
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta_others
  }

  names(beta)[i_gc] <- "gf"
  beta[order(names(beta))]
}
