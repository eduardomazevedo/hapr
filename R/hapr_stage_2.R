#' HAPR second stage fit
#' #'
#' @description
#' Fits the full HAPR model given the first stage fit and an improvement ratio or r2_future.
#'
#' @param first_stage A hapr_first_stage_fit object
#' @param improvement_ratio Ratio to extrapolate by (optional if r2_future supplied)
#' @param r2_current Optional R² of the current fit
#' @param r2_future Optional R² of the future fit (implies improvement_ratio)
#' @return A hapr_fit object
#' Fits the full HAPR model given the first stage fit and an improvement ratio.
#' @export
hapr_second_stage <- function(first_stage,
                              improvement_ratio = NULL,
                              r2_current = NULL,
                              r2_future = NULL) {
  if (!inherits(first_stage, "hapr_first_stage_fit")) {
    stop("first_stage must be a hapr_first_stage_fit object.")
  }

  if (is.null(improvement_ratio) && is.null(r2_future)) {
    stop("Need improvement_ratio or r2_future must be specified.")
  }
  if (!is.null(improvement_ratio) && !is.null(r2_future)) {
    stop("Only one of improvement_ratio or r2_future should be provided.")
  }

  # R² bookkeeping
  if (is.null(r2_current)) {
    r2_current <- tryCatch(as.numeric(first_stage$regressions$y_on_gc$r2), error = function(e) NA_real_)
    if (first_stage$model_type == "cox") r2_current <- NA_real_
  }
  if (is.null(r2_future)) {
    r2_future <- if (is.na(r2_current)) NA_real_ else improvement_ratio * r2_current
  } else {
    improvement_ratio <- r2_future / r2_current
  }

  # tolerate tiny round-off at the boundary
  tol <- 1e-12
  max_ir <- first_stage$stats$max_improvement_ratio
  if (is.finite(max_ir) && improvement_ratio >= max_ir - tol) {
    stop("Improvement ratio too large.")
  }

  # posterior pieces
  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- first_stage$stats$var_v_plus_var_epsilon - var_epsilon
  posterior <- abc(var_epsilon, var_v)

  # map to beta
  beta <- calculate_beta(first_stage$model_type, first_stage$coefficients, posterior)

  # --- Delta method variance ---
  gamma_hat <- first_stage$coefficients$gamma
  theta_hat <- first_stage$coefficients$theta

  Vg <- first_stage$coefficients$vcov_gamma
  Vt <- first_stage$coefficients$vcov_theta
  ng <- length(gamma_hat); nt <- length(theta_hat)

  if (first_stage$model_type %in% c("lm", "probit")) {
    var_v_plus_var_epsilon_hat <- first_stage$stats$var_v_plus_var_epsilon
    vcov_var_v_plus_var_epsilon <- first_stage$regressions$gc_on_w$var_sigma_squared

    param_names <- c(names(gamma_hat), names(theta_hat), "var_v_plus_var_epsilon")
    vcov_full <- matrix(0, ng + nt + 1, ng + nt + 1)
    rownames(vcov_full) <- colnames(vcov_full) <- param_names
    vcov_full[names(gamma_hat), names(gamma_hat)] <- Vg
    vcov_full[names(theta_hat), names(theta_hat)] <- Vt
    vcov_full["var_v_plus_var_epsilon", "var_v_plus_var_epsilon"] <- vcov_var_v_plus_var_epsilon

    J <- calculate_analytical_jacobian(
      model_type = first_stage$model_type,
      gamma = gamma_hat,
      theta = theta_hat,
      var_total = var_v_plus_var_epsilon_hat,
      posterior = posterior,
      beta = beta,
      derived_vars = list(var_epsilon = var_epsilon)
    )
  } else {
    # canonical order for delta method
    param_names <- c(names(gamma_hat), names(theta_hat))
    param_hat   <- stats::setNames(c(gamma_hat, theta_hat), param_names)

    vcov_full <- matrix(0, ng + nt, ng + nt)
    rownames(vcov_full) <- colnames(vcov_full) <- param_names
    vcov_full[names(gamma_hat), names(gamma_hat)] <- Vg
    vcov_full[names(theta_hat), names(theta_hat)] <- Vt

    beta_wrapper <- function(par) {
      gamma <- par[seq_len(ng)]; names(gamma) <- names(gamma_hat)
      theta <- par[ng + seq_len(nt)]; names(theta) <- names(theta_hat)
      calculate_beta(first_stage$model_type, list(gamma = gamma, theta = theta), posterior)
    }

    J <- numDeriv::jacobian(beta_wrapper, param_hat)
    rownames(J) <- names(beta_wrapper(param_hat))
    colnames(J) <- param_names
  }

  vcov_beta <- J %*% vcov_full %*% t(J)
  sd_beta <- sqrt(diag(vcov_beta)); names(sd_beta) <- rownames(J)

  z <- stats::qnorm(0.975)
  ci_beta <- data.frame(
    Estimate  = beta[rownames(J)],
    Std.Error = sd_beta,
    Lower     = beta[rownames(J)] - z * sd_beta,
    Upper     = beta[rownames(J)] + z * sd_beta,
    row.names = rownames(J),
    check.names = FALSE
  )

  # update stripped models (if present)
  if (!is.null(first_stage$regressions$y_on_gf_w$stripped_model)) {
    first_stage$regressions$y_on_gf_w$stripped_model$coefficients <- beta
  }
  first_stage$regressions$y_on_gc_w$coefficients <- beta

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

# ---- internal helper: calculate_beta -----------------------------------------

#' Calculate beta coefficients based on model type
#'
#' Internal helper used by hapr_second_stage().
#' @keywords internal
calculate_beta <- function(model_type, coefficients, posterior) {
  gamma <- coefficients$gamma
  theta <- coefficients$theta

  i_gc    <- which(names(gamma) == "gc")
  i_other <- setdiff(seq_along(gamma), i_gc)

  beta <- gamma

  if (model_type == "lm") {
    beta[i_gc] <- gamma[i_gc] / posterior$a
    theta_others <- theta[names(theta) != "(Intercept)"]
    common <- intersect(names(gamma)[i_other], names(theta_others))
    if (length(common)) {
      beta[i_other][common] <- gamma[i_other][common] -
        posterior$b * theta_others[common] * as.numeric(beta[i_gc])
      if (length(setdiff(names(gamma)[i_other], common))) {
        beta[i_other][setdiff(names(gamma)[i_other], common)] <-
          gamma[i_other][setdiff(names(gamma)[i_other], common)]
      }
    } else {
      beta[i_other] <- gamma[i_other]
    }

  } else if (model_type == "probit") {
    sqrt_input <- posterior$a^2 - (as.numeric(gamma[i_gc])^2) * (posterior$c^2)
    if (sqrt_input < 0) stop("Invalid posterior parameters: sqrt_input is negative")
    beta[i_gc] <- as.numeric(gamma[i_gc]) / sqrt(sqrt_input)
    scale <- sqrt(1 + (posterior$c^2) * (as.numeric(beta[i_gc])^2))
    beta[i_other] <- gamma[i_other] * scale
    theta_others <- theta[names(theta) != "(Intercept)"]
    common <- intersect(names(gamma)[i_other], names(theta_others))
    if (length(common)) {
      beta[i_other][common] <- beta[i_other][common] -
        posterior$b * theta_others[common] * as.numeric(beta[i_gc])
    }

  } else if (model_type == "cox") {
    beta[i_gc] <- gamma[i_gc] / posterior$a
    theta_others <- theta[-1]
    common <- intersect(names(gamma)[i_other], names(theta_others))
    if (length(common)) {
      beta[i_other][common] <- gamma[i_other][common] -
        posterior$b * theta_others[common] * as.numeric(beta[i_gc])
      if (length(setdiff(names(gamma)[i_other], common))) {
        beta[i_other][setdiff(names(gamma)[i_other], common)] <-
          gamma[i_other][setdiff(names(gamma)[i_other], common)]
      }
    } else {
      beta[i_other] <- gamma[i_other]
    }
  }

  names(beta)[i_gc] <- "gf"
  beta[order(names(beta))]
}

#' Calculate analytical Jacobian for delta-method SEs
#'
#' @keywords internal
calculate_analytical_jacobian <- function(model_type, gamma, theta, var_total, posterior, beta, derived_vars) {
  if (!model_type %in% c("lm", "probit")) {
    stop("Analytical Jacobian only implemented for lm and probit.")
  }

  k <- derived_vars$var_epsilon
  x <- var_total
  a <- posterior$a
  b <- posterior$b
  c2 <- posterior$c^2

  da_dx <- k / x^2
  db_dx <- -k / x^2
  dc2_dx <- k^2 / x^2

  gamma_names <- names(gamma)
  theta_names <- names(theta)
  beta_names <- names(beta)

  col_names <- c(gamma_names, theta_names, "var_v_plus_var_epsilon")
  J <- matrix(0, nrow = length(beta_names), ncol = length(col_names))
  rownames(J) <- beta_names
  colnames(J) <- col_names

  if (!"gc" %in% gamma_names) {
    stop("Gamma must include 'gc' for analytical Jacobian.")
  }
  if (!"gf" %in% beta_names) {
    stop("Beta must include 'gf' for analytical Jacobian.")
  }

  row_gf <- which(beta_names == "gf")
  col_gc <- which(gamma_names == "gc")
  col_vt <- length(col_names)

  gamma_gc <- as.numeric(gamma["gc"])
  beta_gc <- as.numeric(beta["gf"])

  gamma_other_names <- setdiff(gamma_names, "gc")
  theta_other_names <- setdiff(theta_names, "(Intercept)")
  common_names <- intersect(gamma_other_names, theta_other_names)

  if (model_type == "lm") {
    dBg_dgg <- 1 / a
    dBg_dx <- -gamma_gc / (a^2) * da_dx

    J[row_gf, col_gc] <- dBg_dgg
    J[row_gf, col_vt] <- dBg_dx

    if (length(gamma_other_names)) {
      for (gname in gamma_other_names) {
        row_g <- which(beta_names == gname)
        col_g <- which(gamma_names == gname)

        J[row_g, col_g] <- 1

        if (gname %in% common_names) {
          theta_g <- as.numeric(theta[gname])
          col_t <- length(gamma_names) + which(theta_names == gname)

          J[row_g, col_gc] <- -theta_g * (b / a)
          J[row_g, col_t] <- -beta_gc * b
          J[row_g, col_vt] <- -theta_g * (db_dx * beta_gc + b * dBg_dx)
        }
      }
    }
  } else if (model_type == "probit") {
    D_sq <- a^2 - (gamma_gc^2) * c2
    if (D_sq <= 0) {
      stop("Invalid posterior parameters: sqrt_input is non-positive")
    }
    D <- sqrt(D_sq)

    dBg_dgg <- (a^2) / (D^3)
    dD_dx <- (1 / D) * (a * da_dx - 0.5 * gamma_gc^2 * dc2_dx)
    dBg_dx <- (-gamma_gc / D_sq) * dD_dx

    J[row_gf, col_gc] <- dBg_dgg
    J[row_gf, col_vt] <- dBg_dx

    S <- sqrt(1 + c2 * beta_gc^2)
    dS_dgg <- (c2 * beta_gc / S) * dBg_dgg
    dS_dx <- (1 / (2 * S)) * (dc2_dx * beta_gc^2 + 2 * c2 * beta_gc * dBg_dx)

    if (length(gamma_other_names)) {
      for (gname in gamma_other_names) {
        row_g <- which(beta_names == gname)
        col_g <- which(gamma_names == gname)
        gamma_g <- as.numeric(gamma[gname])

        J[row_g, col_g] <- S
        J[row_g, col_gc] <- gamma_g * dS_dgg
        J[row_g, col_vt] <- gamma_g * dS_dx

        if (gname %in% common_names) {
          theta_g <- as.numeric(theta[gname])
          col_t <- length(gamma_names) + which(theta_names == gname)

          J[row_g, col_gc] <- J[row_g, col_gc] - b * theta_g * dBg_dgg
          J[row_g, col_t] <- -b * beta_gc
          J[row_g, col_vt] <- J[row_g, col_vt] -
            theta_g * (db_dx * beta_gc + b * dBg_dx)
        }
      }
    }
  }

  J
}
