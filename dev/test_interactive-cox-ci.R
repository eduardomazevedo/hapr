# Create a small interactive test for hapr_cox
set.seed(123)
devtools::load_all()
library(tidyverse)
library(listviewer)
library(survival)

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
  
  linpred <- beta_gf * gf + beta_w1 * w$w1
  baseline_time <- rexp(n, rate = exp(linpred))
  censoring <- rexp(n, rate = 1)
  time <- pmin(baseline_time, censoring)
  event <- as.integer(baseline_time <= censoring)
  
  y <- Surv(time, event)
  
  list(
    y = y,
    gc = gc_normalized,
    w = w,
    true_improvement_ratio = true_improvement_ratio
  )
}

# ---- Part 2: Fit once to get "true" betas ----
sim_data <- create_simulated_dataset()
fit <- hapr(sim_data$y, sim_data$gc, sim_data$w,
            model_type = "cox",
            improvement_ratio = sim_data$true_improvement_ratio)

true_beta <- fit$coefficients$beta
beta_names <- names(true_beta)

# ---- Part 3: Simulations with reused true_beta ----
n_sim <- 1000
covered_matrix <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
fitted_beta_mat <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
colnames(covered_matrix) <- colnames(fitted_beta_mat) <- beta_names

for (i in 1:n_sim) {
  if (i %% 50 == 0) cat("Simulation", i, "\n")
  
  sim_data_i <- create_simulated_dataset(params = list(seed = i))
  fit_i <- hapr(
    y = sim_data_i$y,
    gc = sim_data_i$gc,
    w = sim_data_i$w,
    model_type = "cox",
    improvement_ratio = sim_data_i$true_improvement_ratio
  )
  
  beta_hat_i <- fit_i$coefficients$beta
  ci_beta_i <- fit_i$ci_beta
  
  for (term in beta_names) {
    ci_lower <- ci_beta_i[term, "Lower"]
    ci_upper <- ci_beta_i[term, "Upper"]
    true_val <- true_beta[term]
    
    covered_matrix[i, term] <- (true_val >= ci_lower) & (true_val <= ci_upper)
    fitted_beta_mat[i, term] <- beta_hat_i[term]
  }
}

# ---- Part 4: Summary and histogram plots ----
coverage_df <- tibble(
  Term = beta_names,
  Coverage = round(colMeans(covered_matrix, na.rm = TRUE), 3)
)
print(coverage_df)

for (term in beta_names) {
  estimates <- fitted_beta_mat[, term]
  ci_lower <- as.numeric(fit$ci_beta[term, "Lower"])
  ci_upper <- as.numeric(fit$ci_beta[term, "Upper"])
  true_val <- true_beta[term]
  
  buffer <- 0.05 * (max(estimates) - min(estimates))
  x_min <- min(min(estimates), ci_lower) - buffer
  x_max <- max(max(estimates), ci_upper) + buffer
  
  hist(estimates,
       breaks = 30,
       main = paste("Sampling Distribution of Beta:", term),
       xlab = paste("Estimated beta:", term),
       col = "lightblue",
       border = "white",
       xlim = c(x_min, x_max))
  
  abline(v = ci_lower, col = "red", lwd = 2)
  abline(v = ci_upper, col = "red", lwd = 2)
  abline(v = true_val, col = "blue", lwd = 2, lty = 2)
  
  legend("topright",
         legend = c("CI Lower (delta)", "CI Upper (delta)", "True beta"),
         col = c("red", "red", "blue"),
         lty = c(1, 1, 2),
         lwd = 2,
         bty = "n")
  
  if (interactive()) readline(prompt = "Press [Enter] to continue...")
}