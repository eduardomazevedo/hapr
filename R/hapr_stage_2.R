#' HAPR second stage fit
#'
#' @description
#' Fits the full HARP model given the first stage fit and an improvement ratio.
#' Uses manual Delta Method for standard errors.
#' @param first_stage A hapr_first_stage_fit object
#' @param improvement_ratio The ratio to extrapolate by (required)
#' @return A hapr_lm_fit object containing the results of the second stage
#' @export
hapr_second_stage <- function(
    first_stage,
    improvement_ratio) {

  if (!inherits(first_stage, "hapr_first_stage_fit")) {
    stop("first_stage must be a hapr_first_stage_fit object.")
  }

  if (missing(improvement_ratio) || is.null(improvement_ratio)) {
    stop("improvement_ratio must be specified.")
  }

  # Calculate R-squared values from first stage
  r2_current <- first_stage$regressions$y_on_gc$r2 |> as.numeric()
  r2_future <- improvement_ratio * r2_current

  if (improvement_ratio >= first_stage$stats$max_improvement_ratio) {
    stop(sprintf("Improvement ratio must be less than %s.",
                 first_stage$stats$max_improvement_ratio))
  }

  # --- 1. Point Estimates ---
  derived_parameters <- calculate_parameters(first_stage$model_type, first_stage$parameters, improvement_ratio)

  beta <- derived_parameters$beta
  var_v <- derived_parameters$var_v
  var_epsilon <- derived_parameters$var_epsilon
  posterior <- derived_parameters$posterior

  # --- 2. Manual Delta Method for SEs ---
  
  # Extract inputs
  gamma_hat <- first_stage$parameters$gamma
  theta_hat <- first_stage$parameters$theta
  var_v_plus_var_epsilon_hat <- first_stage$parameters$var_v_plus_var_epsilon
  
  vcov_gamma <- first_stage$vcov_parameters$gamma
  vcov_theta <- first_stage$vcov_parameters$theta
  vcov_var_v_plus_var_epsilon <- first_stage$vcov_parameters$var_v_plus_var_epsilon

  # Construct full covariance matrix of inputs
  # Order: [Gamma, Theta, Var_Total]
  ng <- length(gamma_hat)
  nt <- length(theta_hat)
  
  vcov_full <- matrix(0, ng + nt + 1, ng + nt + 1)
  vcov_full[1:ng, 1:ng] <- vcov_gamma
  vcov_full[(ng + 1):(ng + nt), (ng + 1):(ng + nt)] <- vcov_theta
  vcov_full[(ng + nt + 1), (ng + nt + 1)] <- vcov_var_v_plus_var_epsilon

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
  y_on_gf_w$parameters <- beta
  y_on_gf_w$vcov_parameters <- vcov_beta
  
  # Initialize parameters list for output
  parameters <- list(beta = beta)
  
  # Add model-specific additional parameters to the same list
  if (first_stage$model_type == "lm") {
    # var_eta is estimated as sigma2_y - beta_g^2 * c^2
    # If users want SE for var_eta, we'd need to expand the Jacobian row for it. 
    # For now, we calculate point estimate as requested.
    parameters$var_eta <- first_stage$regressions$y_on_gc_w$sigma_squared - beta["gf"]^2 * posterior$c^2
  }

  result <- list(
    model_type = first_stage$model_type,
    regressions = c(
      first_stage$regressions,
      list(y_on_gf_w = y_on_gf_w)),
    parameters = parameters,
    vcov_parameters = c(
      first_stage$vcov_parameters,
      list(beta = vcov_beta)),
    standard_errors = sd_beta,
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


#' Calculate Analytical Jacobian for HAPR Stage 2
#'
#' @param model_type "lm" or "probit"
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
  
  # Identify indices using strict ordering
  # Note: gamma ordering is: gc, (Intercept), w1, w2, ...
  #      theta ordering is: (Intercept), w1, w2, ...
  # This ordering is guaranteed by hapr_first_stage and is critical for correct alignment
  i_gc <- 1L  # First element of gamma is always gc
  i_w <- 2L:ng  # Remaining elements are (Intercept), w1, w2, ...
  
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
  col_theta    <- (ng + 1):(ng + nt)  # Theta columns start after gamma columns
  col_vt       <- ng + nt + 1
  
  gamma_gc <- gamma[i_gc]
  gamma_w  <- gamma[i_w]  # Includes intercept and w coefficients, in order: (Intercept), w1, w2, ...
  beta_gc  <- beta[i_gc]
  
  if (model_type == "lm") {
    # --- Linear Model Jacobian ---
    # Beta_gc = Gamma_gc / a
    # Beta_w = Gamma_w - Beta_gc * b * theta
    # 
    # Note: theta ordering is (Intercept), w1, w2, ...
    #       gamma_w ordering is (Intercept), w1, w2, ...
    #       These align element-wise: theta[1] aligns with gamma_w[1] (both intercept),
    #       theta[2] aligns with gamma_w[2] (both w1), etc.
    
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
    # Element-wise operation relying on strict ordering:
    # theta[1] (intercept) * b/a, theta[2] (w1) * b/a, etc.
    # This works because theta and gamma_w have matching positional ordering
    J[row_w, col_gamma_gc] <- -theta * (b / a)
    
    # d(B_w)/d(Theta) = -Beta_gc * b * I
    # Sets diagonal relying on strict ordering:
    # row_w[1] (intercept) -> col_theta[1] (theta intercept)
    # row_w[2] (w1) -> col_theta[2] (theta w1), etc.
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
  }
  
  return(J)
}


#' Calculate beta coefficients based on model type
#'
#' @param model_type The type of model ("lm" or "probit")
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

  # Use strict ordering: gamma is gc, (Intercept), w1, w2, ...
  i_gc <- 1L  # First element of gamma is always gc
  i_other <- 2L:length(gamma)  # Remaining elements are (Intercept), w1, w2, ...

  if (model_type == "lm") {
    beta[i_gc] <- gamma[i_gc] / posterior$a
    # Element-wise operation relying on strict ordering:
    # theta ordering: (Intercept), w1, w2, ...
    # gamma[i_other] ordering: (Intercept), w1, w2, ...
    # These align by position: theta[1] with gamma[i_other][1] (both intercept),
    #                          theta[2] with gamma[i_other][2] (both w1), etc.
    beta[i_other] <- gamma[i_other] - beta[i_gc] * posterior$b * theta

  } else if (model_type == "probit") {
    sqrt_input <- posterior$a^2 - (gamma[i_gc]^2) * (posterior$c^2)
    if (sqrt_input < 0) stop("Invalid posterior parameters: sqrt_input is negative")
    beta[i_gc] <- gamma[i_gc] / sqrt(sqrt_input)
    # Element-wise operation relying on strict ordering:
    # theta and gamma[i_other] align by position as (Intercept), w1, w2, ...
    beta[i_other] <- gamma[i_other] * sqrt(1 + (posterior$c^2) * beta[i_gc]^2) -
      posterior$b * theta * beta[i_gc]
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