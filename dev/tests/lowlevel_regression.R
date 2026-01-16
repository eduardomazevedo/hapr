#' Test script for low-level regression functions
#' 
#' Compares fit_lm and fit_probit with standard lm and glm functions
#' Prints comparison tables for coefficients, vcov, and other statistics

# Source the functions (run from project root)
source("R/lowlevel_regression.R")

# Set seed for reproducibility
set.seed(123)

cat(paste(rep("=", 72), collapse = ""), "\n")
cat("Testing Low-Level Regression Functions\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

# ============================================================================
# LINEAR REGRESSION TEST
# ============================================================================

cat("LINEAR REGRESSION TEST\n")
cat(paste(rep("-", 72), collapse = ""), "\n\n")

# Simulate data for linear regression
n <- 100
p <- 3
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- c("x1", "x2", "x3")
beta_true <- c(2, -1, 0.5, 1.5)  # intercept + 3 predictors
y <- beta_true[1] + X %*% beta_true[-1] + rnorm(n, sd = 0.5)

# Fit with low-level function (include intercept in X)
X_with_int <- cbind(`(Intercept)` = 1, X)
fit_low <- fit_lm(y = y, X = X_with_int)

# Fit with standard lm
df_lm <- data.frame(y = y, X)
fit_std <- lm(y ~ ., data = df_lm)

# Compare coefficients
cat("COEFFICIENTS COMPARISON:\n")
coef_comparison <- data.frame(
  Coefficient = names(fit_low$coefficients),
  LowLevel = fit_low$coefficients,
  Standard = coef(fit_std),
  Difference = fit_low$coefficients - coef(fit_std)
)
print(coef_comparison, row.names = FALSE)
cat("\n")

# Compare vcov
cat("VCOV DIAGONAL (VARIANCES) COMPARISON:\n")
vcov_diag <- data.frame(
  Coefficient = names(fit_low$coefficients),
  LowLevel = diag(fit_low$vcov_coefficients),
  Standard = diag(vcov(fit_std)),
  Difference = diag(fit_low$vcov_coefficients) - diag(vcov(fit_std))
)
print(vcov_diag, row.names = FALSE)
cat("\n")

# Compare sigma squared
cat("SIGMA SQUARED COMPARISON:\n")
sigma_comparison <- data.frame(
  Statistic = c("sigma_squared", "var_sigma_squared"),
  LowLevel = c(fit_low$sigma_squared, fit_low$var_sigma_squared),
  Standard = c(summary(fit_std)$sigma^2, NA),
  Difference = c(fit_low$sigma_squared - summary(fit_std)$sigma^2, NA)
)
print(sigma_comparison, row.names = FALSE)
cat("\n")

# Compare R-squared
cat("R-SQUARED COMPARISON:\n")
r2_comparison <- data.frame(
  Method = c("LowLevel", "Standard"),
  R2 = c(fit_low$r2, summary(fit_std)$r.squared),
  Difference = c(0, fit_low$r2 - summary(fit_std)$r.squared)
)
print(r2_comparison, row.names = FALSE)
cat("\n")

# Full VCOV comparison (print first few rows/cols)
cat("VCOV MATRIX COMPARISON (first 3x3):\n")
cat("LowLevel VCOV:\n")
print(fit_low$vcov_coefficients[1:min(3, nrow(fit_low$vcov_coefficients)), 
                                 1:min(3, ncol(fit_low$vcov_coefficients))])
cat("\nStandard VCOV:\n")
print(vcov(fit_std)[1:min(3, nrow(vcov(fit_std))), 
                    1:min(3, ncol(vcov(fit_std)))])
cat("\nDifference:\n")
print(fit_low$vcov_coefficients[1:min(3, nrow(fit_low$vcov_coefficients)), 
                                 1:min(3, ncol(fit_low$vcov_coefficients))] - 
      vcov(fit_std)[1:min(3, nrow(vcov(fit_std))), 
                    1:min(3, ncol(vcov(fit_std)))])
cat("\n\n")

# ============================================================================
# PROBIT REGRESSION TEST
# ============================================================================

cat("PROBIT REGRESSION TEST\n")
cat(paste(rep("-", 72), collapse = ""), "\n\n")

# Simulate data for probit regression
n_probit <- 200
p_probit <- 2
X_probit <- matrix(rnorm(n_probit * p_probit), nrow = n_probit, ncol = p_probit)
colnames(X_probit) <- c("x1", "x2")
beta_probit_true <- c(-0.5, 1, -0.8)  # intercept + 2 predictors
linear_pred <- beta_probit_true[1] + X_probit %*% beta_probit_true[-1]
prob <- pnorm(linear_pred)
y_probit <- rbinom(n_probit, size = 1, prob = prob)

# Fit with low-level function (include intercept in X)
X_probit_with_int <- cbind(`(Intercept)` = 1, X_probit)
fit_probit_low <- fit_probit(y = y_probit, X = X_probit_with_int)

# Fit with standard glm
df_probit <- data.frame(y = y_probit, X_probit)
fit_probit_std <- glm(y ~ ., data = df_probit, family = binomial(link = "probit"))

# Compare coefficients
cat("COEFFICIENTS COMPARISON:\n")
coef_probit_comparison <- data.frame(
  Coefficient = names(fit_probit_low$coefficients),
  LowLevel = fit_probit_low$coefficients,
  Standard = coef(fit_probit_std),
  Difference = fit_probit_low$coefficients - coef(fit_probit_std)
)
print(coef_probit_comparison, row.names = FALSE)
cat("\n")

# Compare vcov
cat("VCOV DIAGONAL (VARIANCES) COMPARISON:\n")
vcov_probit_diag <- data.frame(
  Coefficient = names(fit_probit_low$coefficients),
  LowLevel = diag(fit_probit_low$vcov_coefficients),
  Standard = diag(vcov(fit_probit_std)),
  Difference = diag(fit_probit_low$vcov_coefficients) - diag(vcov(fit_probit_std))
)
print(vcov_probit_diag, row.names = FALSE)
cat("\n")

# Compare R-squared
cat("R-SQUARED COMPARISON:\n")
cat("Note: LowLevel uses liability-scale R2, Standard uses McFadden's pseudo-R2\n")
# Calculate McFadden's pseudo-R2 for comparison
null_model <- glm(y ~ 1, data = df_probit, family = binomial(link = "probit"))
mcfadden_r2 <- 1 - (logLik(fit_probit_std) / logLik(null_model))
r2_probit_comparison <- data.frame(
  Method = c("LowLevel (liability-scale)", "Standard (McFadden)"),
  R2 = c(fit_probit_low$r2, as.numeric(mcfadden_r2)),
  Note = c("Var(Xb)/(Var(Xb)+1)", "1 - logLik(model)/logLik(null)")
)
print(r2_probit_comparison, row.names = FALSE)
cat("\n")

# Full VCOV comparison
cat("VCOV MATRIX COMPARISON:\n")
cat("LowLevel VCOV:\n")
print(fit_probit_low$vcov_coefficients)
cat("\nStandard VCOV:\n")
print(vcov(fit_probit_std))
cat("\nDifference:\n")
print(fit_probit_low$vcov_coefficients - vcov(fit_probit_std))
cat("\n\n")

# ============================================================================
# TEST WITHOUT INTERCEPT
# ============================================================================

cat("LINEAR REGRESSION TEST (NO INTERCEPT)\n")
cat(paste(rep("-", 72), collapse = ""), "\n\n")

# Fit without intercept
fit_low_no_int <- fit_lm(y = y, X = X)
fit_std_no_int <- lm(y ~ 0 + ., data = df_lm)

cat("COEFFICIENTS COMPARISON:\n")
coef_no_int <- data.frame(
  Coefficient = names(fit_low_no_int$coefficients),
  LowLevel = fit_low_no_int$coefficients,
  Standard = coef(fit_std_no_int),
  Difference = fit_low_no_int$coefficients - coef(fit_std_no_int)
)
print(coef_no_int, row.names = FALSE)
cat("\n")

cat("VCOV DIAGONAL COMPARISON:\n")
vcov_no_int_diag <- data.frame(
  Coefficient = names(fit_low_no_int$coefficients),
  LowLevel = diag(fit_low_no_int$vcov_coefficients),
  Standard = diag(vcov(fit_std_no_int)),
  Difference = diag(fit_low_no_int$vcov_coefficients) - diag(vcov(fit_std_no_int))
)
print(vcov_no_int_diag, row.names = FALSE)
cat("\n")

cat(paste(rep("=", 72), collapse = ""), "\n")
cat("Testing Complete\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
