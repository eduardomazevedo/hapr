# Create a small interactive test for hapr_lm
set.seed(123)
devtools::load_all()
library(tidyverse)
library(listviewer)

# ---- Part 1: Data generation ----
default_params <- list(
  n = 1e4,
  var_v = 1 / 3,
  var_epsilon = 1 / 3,
  var_thetaw = 1 / 3,
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
n_sim <- 1000
covered_matrix <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
fitted_beta_mat <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
colnames(covered_matrix) <- colnames(fitted_beta_mat) <- beta_names

for (i in 1:n_sim) {
  if (i %% 200 == 0) cat("Simulation", i, "\n")
  sim_data <- create_simulated_dataset(params = list(seed = i))

  fit <- hapr(
    y = sim_data$y,
    gc = sim_data$gc_normalized,
    w = sim_data$w,
    model_type = "lm",
    improvement_ratio = sim_data$true_improvement_ratio
  )

  beta_hat <- fit$coefficients$beta
  ci_beta <- fit$ci_beta

  for (term in beta_names) {
    ci_lower <- ci_beta[term, "Lower"]
    ci_upper <- ci_beta[term, "Upper"]
    covered_matrix[i, term] <- (true_beta[term] >= ci_lower) & (true_beta[term] <= ci_upper)
    fitted_beta_mat[i, term] <- beta_hat[term]
  }
}

# ---- Part 4: Coverage and histogram ----
coverage_df <- tibble(
  Term = beta_names,
  Coverage = round(colMeans(covered_matrix, na.rm = TRUE), 3)
)
print(coverage_df)

for (term in beta_names) {
  estimates <- fitted_beta_mat[, term]
  ci_bounds <- quantile(estimates, probs = c(0.05, 0.95), na.rm = TRUE)
  true_val <- true_beta[term]

  buffer <- 0.05 * (max(estimates) - min(estimates))
  x_min <- min(ci_bounds[1], true_val, min(estimates)) - buffer
  x_max <- max(ci_bounds[2], true_val, max(estimates)) + buffer

  hist(estimates,
       breaks = 30,
       main = paste("Sampling Distribution of Beta:", term),
       xlab = paste("Estimated beta:", term),
       col = "lightblue",
       border = "white",
       xlim = c(x_min, x_max))

  abline(v = ci_bounds[1], col = "red", lwd = 2)
  abline(v = ci_bounds[2], col = "red", lwd = 2)
  abline(v = true_val, col = "blue", lwd = 2, lty = 2)

  legend("topright",
         legend = c("CI Lower", "CI Upper", "True beta"),
         col = c("red", "red", "blue"),
         lty = c(1, 1, 2),
         lwd = 2,
         bty = "n")

  if (interactive()) readline(prompt = "Press [Enter] to continue...")
}
