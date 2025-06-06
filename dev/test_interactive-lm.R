# Create a small interactive test for hapr_lm
set.seed(123) # For reproducibility
devtools::load_all()
library(tidyverse)
library(listviewer)

# Create fake data
n <- 1e4

var_v <- 1 / 3
var_epsilon <- 1 / 3
var_thetaw <- 1 / 3

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

y <- 0.42 * gf + rnorm(n) + 0.17 * w$w1


# Call the hapr_lm function
fit <- hapr(y, gc_normalized, w, model_type = "lm", improvement_ratio = true_improvement_ratio)

print(fit$coefficients$beta)

simulated_w <- hapr_simulate(fit, w = w) |> as_tibble()
simulated_w_gc <- hapr_simulate(fit, w = w, gc = gc_normalized) |> as_tibble()

# Test predict
predicted_w <- predict.hapr_fit(fit, newdata = simulated_w)

plot(predicted_w$y_hat_w, predicted_w$y_hat_w)
plot(predicted_w$y_hat_w, predicted_w$y_hat_gc_w)
plot(predicted_w$y_hat_w, predicted_w$y_hat_gf_w)

print.hapr_fit(fit)


# ---- CI coverage test for beta ----

n_sim <- 1000
true_beta <- fit$coefficients$beta
beta_names <- names(true_beta)
covered_matrix <- matrix(NA, nrow = n_sim, ncol = length(true_beta))
colnames(covered_matrix) <- beta_names

check_coverage <- function(fit_sim, beta_true) {
  ci <- fit_sim$ci_beta
  ci <- ci[beta_names, , drop = FALSE]
  (beta_true >= ci$Lower) & (beta_true <= ci$Upper)
}

for (i in 1:n_sim) {
  print(i)
  
  # simulate new gf, gc, w
  sim_data <- hapr_simulate(fit, w = w)
  
  # rebuild y using same model logic
  y_sim <- 0.42 * sim_data$gf + rnorm(n) + 0.17 * sim_data$w1
  gc_sim <- sim_data$gc
  w_sim <- sim_data %>% select(all_of(names(w)))
  
  stopifnot(length(y_sim) == length(gc_sim), length(gc_sim) == nrow(w_sim))
  
  fit_sim <- hapr(y_sim, gc_sim, w_sim, model_type = "lm", improvement_ratio = true_improvement_ratio)
  covered_matrix[i, ] <- check_coverage(fit_sim, true_beta)
}

# Summarize empirical coverage
coverage_df <- colMeans(covered_matrix, na.rm = TRUE) |> round(3)
print(tibble(Term = beta_names, Coverage = coverage_df))

