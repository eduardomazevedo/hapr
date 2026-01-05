#' HAPR second stage fit
#'
#' @description
#' Fits the full HARP model given the first stage fit and an improvement ratio.
#' Uses manual Delta Method for standard errors.
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

  # --- 1. Point Estimates ---
  derived_parameters <- calculate_parameters(first_stage$model_type, first_stage$coefficients, improvement_ratio)

  beta <- derived_parameters$beta
  var_v <- derived_parameters$var_v
  var_epsilon <- derived_parameters$var_epsilon
  posterior <- derived_parameters$posterior

  # --- 2. Manual Delta Method for SEs ---
  
  # Extract inputs
  gamma_hat <- first_stage$coefficients$gamma
  theta_hat <- first_stage$coefficients$theta
  var_v_plus_var_epsilon_hat <- first_stage$coefficients$var_v_plus_var_epsilon
  
  vcov_gamma <- first_stage$vcov_coefficients$gamma
  vcov_theta <- first_stage$vcov_coefficients$theta
  vcov_var_v_plus_var_epsilon <- first_stage$vcov_coefficients$var_v_plus_var_epsilon

  # Validate names
  names_gamma <- names(gamma_hat)
  names_theta <- names(theta_hat)
  stopifnot(all.equal(names_gamma, rownames(vcov_gamma)))
  stopifnot(all.equal(names_theta, rownames(vcov_theta)))

  # Construct full covariance matrix of inputs
  # Order: [Gamma, Theta, Var_Total]
  ng <- length(gamma_hat)
  nt <- length(theta_hat)
  
  # Names for parameter vector
  param_names <- c(names_gamma, names_theta, "var_v_plus_var_epsilon")
  
  vcov_full <- matrix(0, ng + nt + 1, ng + nt + 1)
  vcov_full[1:ng, 1:ng] <- vcov_gamma
  vcov_full[(ng + 1):(ng + nt), (ng + 1):(ng + nt)] <- vcov_theta
  vcov_full[(ng + nt + 1), (ng + nt + 1)] <- vcov_var_v_plus_var_epsilon
  rownames(vcov_full) <- colnames(vcov_full) <- param_names

  # Calculate Analytical Jacobian
  J <- calculate_analytical_jacobian(
    model_type = first_stage$model_type,
    gamma = gamma_hat,
    theta = theta_hat,
    var_total = var_v_plus_var_epsilon_hat,
    posterior = posterior,
    beta = beta,
    derived_vars = list(var_epsilon = var_epsilon)
  )
  
  # Ensure Jacobian names match vcov_full
  # Jacobian columns should correspond to input parameters
  if (!all(colnames(J) == param_names)) {
    # Reorder if necessary (though construction is deterministic)
    J <- J[, param_names, drop = FALSE]
  }

  # Calculate Output Covariance: V_beta = J * V_in * J'
  vcov_beta <- J %*% vcov_full %*% t(J)
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

  # --- 3. Construct Output Objects ---

  # Create model for y_on_gf_w (copy from y_on_gc_w)
  y_on_gf_w <- first_stage$regressions$y_on_gc_w

  # Update model coefficients
  y_on_gf_w$coefficients <- beta
  y_on_gf_w$vcov_coefficients <- vcov_beta
  
  # Update names: change "gc" to "gf"
  names(y_on_gf_w$coefficients)[names(y_on_gf_w$coefficients) == "gc"] <- "gf"
  rownames(y_on_gf_w$vcov_coefficients)[rownames(y_on_gf_w$vcov_coefficients) == "gc"] <- "gf"
  colnames(y_on_gf_w$vcov_coefficients)[colnames(y_on_gf_w$vcov_coefficients) == "gc"] <- "gf"

  additional_parameters <- list()
  if (first_stage$model_type == "lm") {
    # var_eta is estimated as sigma2_y - beta_g^2 * c^2
    # If users want SE for var_eta, we'd need to expand the Jacobian row for it. 
    # For now, we calculate point estimate as requested.
    additional_parameters$var_eta <- first_stage$regressions$y_on_gc_w$sigma_squared - beta["gf"]^2 * posterior$c^2
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


#' Calculate Analytical Jacobian for HAPR Stage 2
#'
#' @param model_type "lm", "probit", or "cox"
#' @param gamma Vector of stage 2 gamma coefficients (includes gc)
#' @param theta Vector of stage 1 theta coefficients
#' @param var_total Scalar, total variance of Gc residuals
#' @param posterior List with a, b, c
#' @param beta Vector of resulting beta coefficients
#' @param derived_vars List containing var_epsilon
#'
#' @return Matrix J of dimensions (length(beta)) x (length(gamma) + length(theta) + 1)
calculate_analytical_jacobian <- function(model_type, gamma, theta, var_total, posterior, beta, derived_vars) {
  
  ng <- length(gamma)
  nt <- length(theta)
  nb <- length(beta)
  
  # Identify indices
  i_gc <- which(names(gamma) == "gc")
  if(length(i_gc) == 0) stop("Could not find 'gc' in gamma coefficients")
  i_w <- which(names(gamma) != "gc")
  
  # Constants/Derivatives w.r.t var_total (x)
  # k = var_epsilon
  # a = 1 - k/x  => da/dx = k/x^2
  x <- var_total
  k <- derived_vars$var_epsilon
  
  da_dx <- k / x^2
  db_dx <- -k / x^2  # b = k/x
  
  # c^2 = k * (1 - k/x) = k - k^2/x
  # dc^2/dx = k^2 / x^2
  dc2_dx <- k^2 / x^2
  
  c2 <- posterior$c^2
  a <- posterior$a
  b <- posterior$b
  
  # Initialize Jacobian
  # Rows: Beta (named same as beta)
  # Cols: Gamma (ng), Theta (nt), VarTotal (1)
  J <- matrix(0, nrow = nb, ncol = ng + nt + 1)
  
  col_names <- c(names(gamma), names(theta), "var_v_plus_var_epsilon")
  colnames(J) <- col_names
  rownames(J) <- names(beta)
  
  # Output indices
  row_gc <- i_gc # Usually 1
  row_w  <- i_w
  
  # Input indices (Columns)
  col_gamma_gc <- i_gc
  col_gamma_w  <- i_w
  col_theta    <- (ng + 1):(ng + nt)
  col_vt       <- ng + nt + 1
  
  gamma_gc <- gamma[i_gc]
  gamma_w  <- gamma[i_w]
  beta_gc  <- beta[i_gc]
  
  if (model_type == "lm") {
    # --- Linear Model Jacobian ---
    # Beta_gc = Gamma_gc / a
    # Beta_w = Gamma_w - Beta_gc * b * theta
    
    # 1. d(Beta_gc) row
    # d(B_gc)/d(G_gc) = 1/a
    J[row_gc, col_gamma_gc] <- 1 / a
    
    # d(B_gc)/d(VarTotal) = -G_gc/a^2 * da/dx
    J[row_gc, col_vt] <- -gamma_gc / (a^2) * da_dx
    
    # 2. d(Beta_w) rows
    # d(B_w)/d(G_w) = Identity
    # Note: Using diagonal assignment logic carefully with names
    J[row_w, col_gamma_w] <- diag(length(row_w))
    
    # d(B_w)/d(G_gc) = -theta * b * d(B_gc)/d(G_gc) = -theta * b/a
    # theta must match dimensions of w.
    # In LM: stopifnot(all(names(theta) == names(gamma[i_other])))
    # So theta maps 1-to-1 with rows of w
    J[row_w, col_gamma_gc] <- -theta * (b / a)
    
    # d(B_w)/d(Theta) = -Beta_gc * b * I
    # This sets the diagonal of the theta-block for B_w rows
    # J[row_w, col_theta] is a square matrix if dims match
    diag_indices <- cbind(row_w, col_theta)
    J[diag_indices] <- -beta_gc * b
    
    # d(B_w)/d(VarTotal)
    # d(B_w)/dx = -theta * [ Beta_gc * db/dx + b * d(B_gc)/dx ]
    #           = -theta * [ Beta_gc * (-k/x^2) + b * (-G_gc/a^2 * k/x^2) ]
    # Simplified in theory files: theta * beta_gc * (k/x^2) * (1/a)
    J[row_w, col_vt] <- theta * beta_gc * (k / x^2) * (1/a)
    
  } else if (model_type == "probit") {
    # --- Probit Jacobian ---
    # Beta_gc = Gamma_gc / D, where D = sqrt(a^2 - Gamma_gc^2 * c^2)
    # Beta_w = Gamma_w * S - b * theta * Beta_gc, where S = sqrt(1 + c^2 * Beta_gc^2)
    
    D <- gamma_gc / beta_gc # recovered D
    D_sq <- D^2
    
    # 1. d(Beta_gc) row
    # d(B_gc)/d(G_gc) = a^2 / D^3
    dBg_dgg <- (a^2) / (D^3)
    J[row_gc, col_gamma_gc] <- dBg_dgg
    
    # d(B_gc)/d(VarTotal)
    # dD/dx = (1/D) * (a * da_dx - 0.5 * G_gc^2 * dc2_dx)
    dD_dx <- (1/D) * (a * da_dx - 0.5 * gamma_gc^2 * dc2_dx)
    dBg_dx <- (-gamma_gc / D_sq) * dD_dx
    J[row_gc, col_vt] <- dBg_dx
    
    # 2. d(Beta_w) rows
    S <- sqrt(1 + c2 * beta_gc^2)
    
    # d(B_w)/d(G_w) = S * I
    J[row_w, col_gamma_w] <- diag(length(row_w)) * S
    
    # d(B_w)/d(Theta) = -b * Beta_gc * I
    diag_indices <- cbind(row_w, col_theta)
    J[diag_indices] <- -b * beta_gc
    
    # d(B_w)/d(G_gc)
    # dS/dG_gc = (c2 * B_gc / S) * d(B_gc)/d(G_gc)
    dS_dgg <- (c2 * beta_gc / S) * dBg_dgg
    # d(B_w)/dG_gc = G_w * dS/dG_gc - b * theta * d(B_gc)/dG_gc
    # Note: This is a vector assignment (length of w)
    J[row_w, col_gamma_gc] <- gamma_w * dS_dgg - (b * theta * dBg_dgg)
    
    # d(B_w)/d(VarTotal)
    # dS/dx = (1/2S) * (dc2_dx * B_gc^2 + 2 * c2 * B_gc * d(B_gc)/dx)
    dS_dx <- (1/(2*S)) * (dc2_dx * beta_gc^2 + 2 * c2 * beta_gc * dBg_dx)
    # d(B_w)/dx = G_w * dS/dx - theta * (db_dx * B_gc + b * d(B_gc)/dx)
    J[row_w, col_vt] <- gamma_w * dS_dx - theta * (db_dx * beta_gc + b * dBg_dx)
    
  } else if (model_type == "cox") {
    # --- Cox Jacobian ---
    # Structure similar to LM, but handling Theta/Gamma misalignment
    # Logic from calculate_parameters:
    # theta_intercept <- theta[1]
    # theta_others <- theta[-1]
    # beta_w = gamma_w - beta_gc * b * theta_others
    
    theta_others <- theta[-1]
    
    # 1. d(Beta_gc) row (Same as LM)
    J[row_gc, col_gamma_gc] <- 1 / a
    J[row_gc, col_vt] <- -gamma_gc / (a^2) * da_dx
    
    # 2. d(Beta_w) rows
    # d(B_w)/d(G_w) = Identity
    J[row_w, col_gamma_w] <- diag(length(row_w))
    
    # d(B_w)/d(G_gc) = -theta_others * b/a
    J[row_w, col_gamma_gc] <- -theta_others * (b / a)
    
    # d(B_w)/d(Theta)
    # Theta has intercept (index 1) and others (indices 2..nt)
    # B_w does not depend on Theta_intercept
    # d(B_w)/d(Theta_others) = -Beta_gc * b * I
    
    # We map row_w to col_theta[2...nt]
    # Check dimensions
    if (length(row_w) != length(theta_others)) {
       stop("Dimension mismatch in Cox Jacobian: Gamma_w vs Theta_others")
    }
    
    # Identify the columns for theta_others within the full J
    col_theta_others <- col_theta[-1]
    diag_indices <- cbind(row_w, col_theta_others)
    J[diag_indices] <- -beta_gc * b
    
    # d(B_w)/d(VarTotal)
    # Same form as LM, but using theta_others
    J[row_w, col_vt] <- theta_others * beta_gc * (k / x^2) * (1/a)
  }
  
  return(J)
}


#' Calculate beta coefficients based on model type
#'
#' @param model_type The type of model ("lm", "probit" or "cox")
#' @param coefficients A list containing gamma and theta coefficients
#' @param improvement_ratio The ratio to extrapolate by
#' @return A list with posterior stats and beta coefficients
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
  names(beta)[i_other] <- names(gamma)[i_other]
  
  list(
    var_epsilon = var_epsilon,
    var_v = var_v,
    posterior = posterior,
    beta = beta)
}

#' Helper to calculate a, b, c posterior parameters
abc <- function(var_epsilon, var_v) {
  var_total <- var_v + var_epsilon
  a <- var_v / var_total
  b <- var_epsilon / var_total
  c_val <- sqrt(var_epsilon * var_v / var_total)
  list(a = a, b = b, c = c_val)
}