# Simulating survival data with a hazard rate dependent on a latent variable.
# gf and gc are normally distributed, and the survival time follows an
# exponential distribution with a hazard rate of exp(gc) / 10.

set.seed(123) # For reproducibility
library(survival)
devtools::load_all()
library(tidyverse)
library(listviewer)

# Number of observations
n <- 1e3 # Change as needed
var_v <- 1 / 3
var_epsilon <- 2 / 3

true_improvement_ratio <- 1 / (1 - var_epsilon)

# Simulating gf ~ N(0, 1/3)
gf <- rnorm(n, mean = 0, sd = sqrt(var_v))

# Simulating w ~ N(0, 1)
w <- rnorm(n, mean = 0, sd = 1) |> as_tibble()
names(w) <- "w1"

# Simulating epsilon ~ N(0, 2/3)
epsilon <- rnorm(n, mean = 0, sd = sqrt(var_epsilon))

# Computing gc = gf + epsilon
gc <- gf + epsilon

# Hazard rate: exp(0.42 * gf + 0.17 * w)
hazard_rate <- exp(0.42 * gf + 0.17 * w$w1)

# Survival time t ~ Exp(hazard_rate)
t <- rexp(n, rate = hazard_rate)

# Storing data in a dataframe
sim_data <- data.frame(gf, w$w1, gc, hazard_rate, t)

# Display the first few rows
head(sim_data)

# Fit
fit <- hapr(Surv(sim_data$t), gc, w, model_type = "cox", improvement_ratio = true_improvement_ratio)

# Print the results
# str(fit)

listviewer::jsonedit(fit)

print(fit$coefficients$beta)

simulated_w <- hapr_simulate(fit, w = w) |> as_tibble()

# Test predict
predicted_w <- predict(fit, newdata = simulated_w)
summary(predicted_w)

# Test basehaz
basehaz <- hapr_basehaz(fit, covariates = "gf_w")
summary(basehaz)

# Test survfit
survival_curves <- hapr_survfit(fit, covariates = "gf_w", newdata = simulated_w |> head(300))

# Test plot.hapr_survfit
plot(survival_curves, mode = "percentiles")


# Conditional survival curves starting at time = 2
survival_cond <- hapr_survfit(
  fit,
  covariates = "gf_w",
  newdata = simulated_w |> head(300),
  start.time = 2
)

plot(survival_cond, mode = "percentiles", percentiles = c(0.1, 0.5, 0.9))



# Predict relative risk
pred_risk <- predict(fit, newdata = simulated_w, covariates = "gf_w", type = "risk")$y_hat_gf_w
threshold <- quantile(pred_risk, 0.9)

# Compute average survival for top 10% risk
sf_avg <- hapr_survfit(
  fit,
  covariates = "gf_w",
  newdata = simulated_w,
  aggregate = pred_risk > threshold
)

plot(sf_avg$time, sf_avg$surv_avg, type = "l", lwd = 2,
     main = "Average Survival (Top 10% Risk)", xlab = "Time", ylab = "Survival")


sf_ci <- hapr_survfit(
  fit,
  covariates = "gf_w",
  newdata = simulated_w,
  aggregate = pred_risk > threshold,
  conf.int = TRUE,
  conf.level = 0.95,
  n.boot = 500
)

# Plot with bands
plot(sf_ci$time, sf_ci$surv_avg, type = "l", ylim = c(0, 1),
     xlab = "Time", ylab = "Survival", main = "Avg Survival + 95% CI (Top 10%)", lwd = 2)
lines(sf_ci$time, sf_ci$surv_avg_lower, col = "red", lty = 2)
lines(sf_ci$time, sf_ci$surv_avg_upper, col = "red", lty = 2)
legend("topright", legend = c("Mean", "95% CI"), col = c("black", "red"), lty = c(1, 2))

