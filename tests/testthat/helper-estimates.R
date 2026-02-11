make_true_coef <- function(beta_g, beta_w) {
  true_coef <- c(beta_g, beta_w)
  names(true_coef) <- c("gf", "(Intercept)", paste0("w", seq_len(length(beta_w) - 1)))
  true_coef
}

align_ci_beta <- function(ci_beta, true_coef) {
  coef_names <- rownames(ci_beta)
  estimates <- setNames(ci_beta$Estimate, coef_names)
  se <- setNames(ci_beta$Std.Error, coef_names)
  lower_ci <- setNames(ci_beta$Lower, coef_names)
  upper_ci <- setNames(ci_beta$Upper, coef_names)

  list(
    estimates = estimates[names(true_coef)],
    se = se[names(true_coef)],
    lower_ci = lower_ci[names(true_coef)],
    upper_ci = upper_ci[names(true_coef)]
  )
}
