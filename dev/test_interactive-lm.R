# Create a small interactive test for hapr_lm
set.seed(123)
devtools::load_all()
library(tidyverse)

# ---- Part 1: Data generation ----
default_params <- list(
  n = 1e4,
  var_v = 1 / 8,
  var_epsilon = 3 / 4,
  var_thetaw = 1 / 8,
  beta_gf = 0.42,
  beta_w1 = 0.17,
  seed = 123
)

create_simulated_dataset <- function(params = list()) {
  p <- modifyList(default_params, params)
  set.seed(p$seed)

  n <- p$n
  var_v <- p$var_v
  var_epsilon <- p$var_epsilon
  var_thetaw <- p$var_thetaw
  beta_gf <- p$beta_gf
  beta_w1 <- p$beta_w1

  true_improvement_ratio <- 1 / (1 - var_epsilon)

  w <- data.frame(
    w1 = rnorm(n),
    w2 = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )

  v <- rnorm(n) * sqrt(var_v)
  epsilon <- rnorm(n) * sqrt(var_epsilon)

  gf <- w$w1 * sqrt(var_thetaw) + v
  gc <- gf + epsilon
  gc_normalized <- scale(gc) |> as.numeric()

  y <- beta_gf * gf + rnorm(n) + beta_w1 * w$w1

  list(
    y = y,
    gc_normalized = gc_normalized,
    w = w,
    true_improvement_ratio = true_improvement_ratio
  )
}

# ---- Part 2: Set true beta explicitly ----
true_beta <- c(gf = 0.42, w1 = 0.17, w2B = 0, w2C = 0)
beta_names <- names(true_beta)

# ---- Part 3: Run simulations ----
n_sim <- 100
covered_matrix <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
fitted_beta_mat <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
se_beta_mat <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
colnames(covered_matrix) <- colnames(fitted_beta_mat) <- colnames(se_beta_mat) <- beta_names

for (i in 1:n_sim) {
  if (i %% 20 == 0) cat("Simulation", i, "\n")
  sim_data <- create_simulated_dataset(params = list(seed = i, n = 10000))

  fit <- hapr(
    y = sim_data$y,
    gc = sim_data$gc_normalized,
    w = sim_data$w,
    model_type = "lm",
    improvement_ratio = sim_data$true_improvement_ratio
  )

  beta_hat <- fit$coefficients$beta
  ci_beta <- fit$ci_beta
  se_hat <- fit$standard_errors

  for (term in beta_names) {
    ci_lower <- ci_beta[term, "Lower"]
    ci_upper <- ci_beta[term, "Upper"]
    covered_matrix[i, term] <- (true_beta[term] >= ci_lower) & (true_beta[term] <= ci_upper)
    fitted_beta_mat[i, term] <- beta_hat[term]
    se_beta_mat[i, term] <- se_hat[term]
  }
}

# ---- Part 4: Summary Table ----
summary_table <- tibble(
  Parameter = beta_names,
  True = true_beta[beta_names],
  Avg_Estimate = colMeans(fitted_beta_mat, na.rm = TRUE),
  Avg_SE = colMeans(se_beta_mat, na.rm = TRUE),
  Stdev_Estimates = apply(fitted_beta_mat, 2, sd, na.rm = TRUE),
  Coverage_Pct = round(100 * colMeans(covered_matrix, na.rm = TRUE), 2)
)

print(summary_table)
