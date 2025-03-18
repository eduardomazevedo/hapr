# Simulating survival data with a hazard rate dependent on a latent variable.
# gf and gc are normally distributed, and the survival time follows an 
# exponential distribution with a hazard rate of exp(gc) / 10.

set.seed(123)  # For reproducibility
library(survival)
devtools::load_all()
library(tidyverse)
library(listviewer)

# Number of observations
n <- 1e3  # Change as needed
var_v <- 1/3
var_epsilon <- 2/3

true_improvement_ratio <- 1 / (1 - var_epsilon)

# Simulating gf ~ N(0, 1/3)
gf <- rnorm(n, mean = 0, sd = sqrt(var_v))

# Simulating w ~ N(0, 1)
w <- rnorm(n, mean = 0, sd = 1) |> as.tibble()
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
fit <- hapr(Surv(sim_data$t), gc, w , model_type = "cox", improvement_ratio = true_improvement_ratio)

# Print the results
# str(fit)

listviewer::jsonedit(fit)

print(fit$coefficients$beta)

simulated_w <- simulate(fit, w = w) |> as.tibble()

# Test predict
predicted_w <- predict.hapr_fit(fit, newdata = simulated_w)
summary(predicted_w)

# Test basehaz
basehaz <- basehaz.hapr_fit(fit, covariates = "gf_w")
summary(basehaz)

# Test survfit
survfit <- survfit.hapr_fit(fit, covariates = "gf_w", newdata = simulated_w |> head(300))

# Test plot.hapr_survfit
plot.hapr_survfit(survfit, mode = "percentiles")
# plot.hapr_survfit(survfit, mode = "subjects", newdata = simulated_w |> head(3))