stage2_params_default <- function() {
  list(
    beta_g = 1.42,
    beta_w = c(0.1, 0.17, 0.27, -0.27),
    theta = c(0.0, 0.1, -0.2, 0.3),
    var_y = 1.0
  )
}

survival_params_default <- function() {
  list(
    beta_g = 0.6,
    beta_w = c(0.1, -0.2, 0.15, 0.05),
    theta = c(0.0, 0.1, -0.25, 0.2),
    censor_rate = 0.2
  )
}
