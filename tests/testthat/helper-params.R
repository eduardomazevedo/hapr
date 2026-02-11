stage2_params_default <- function() {
  list(
    beta_g = 1.42,
    beta_w = c(0.1, -0.2, 0.3, -0.4),
    theta = c(0.0, 1.0, -2.0, 3.0),
    var_y = 1.0
  )
}

survival_params_default <- function() {
  list(
    beta_g = 1.42,
    beta_w = c(0.1, -0.2, 0.3, -0.4),
    theta = c(0.0, 1.0, -2.0, 3.0),
    censor_rate = 0.2
  )
}

normalize_theta <- function(theta, var_v, var_epsilon) {
  target_var <- 1 - var_v - var_epsilon
  if (target_var < 0) {
    stop("target variance for theta'w must be non-negative")
  }
  current_var <- sum(theta[-1]^2)
  if (current_var <= 0) {
    stop("theta must have non-zero variance contribution in non-intercept terms")
  }
  scale <- sqrt(target_var / current_var)
  theta_scaled <- theta
  theta_scaled[-1] <- theta_scaled[-1] * scale
  theta_scaled
}
