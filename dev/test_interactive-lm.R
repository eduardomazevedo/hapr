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
