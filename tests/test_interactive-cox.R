# Simulating survival data with a hazard rate dependent on a latent variable.
# gf and gc are normally distributed, and the survival time follows an 
# exponential distribution with a hazard rate of exp(gc) / 10.

set.seed(123)  # For reproducibility

# Number of observations
n <- 1e4  # Change as needed
var_v <- 1/3
var_epsilon <- 2/3

true_improvement_ratio <- 1 / (1 - var_epsilon)

# Simulating gf ~ N(0, 1/3)
gf <- rnorm(n, mean = 0, sd = sqrt(var_v))

# Simulating w ~ N(0, 1)
w <- rnorm(n, mean = 0, sd = 1)

# Simulating epsilon ~ N(0, 2/3)
epsilon <- rnorm(n, mean = 0, sd = sqrt(var_epsilon))

# Computing gc = gf + epsilon
gc <- gf + epsilon

# Hazard rate: exp(0.42 * gf + 0.17 * w) / 10
hazard_rate <- exp(0.42 * gf + 0.17 * w) / 10

# Survival time t ~ Exp(hazard_rate)
t <- rexp(n, rate = hazard_rate)

# Storing data in a dataframe
sim_data <- data.frame(gf, w, gc, hazard_rate, t)

# Display the first few rows
head(sim_data)

# Now fit the model
fit <- hapr_cox(Surv(sim_data$t), sim_data$gc, sim_data$w |> as.data.frame(), improvement_ratio = true_improvement_ratio)

# Print the results
print(fit)

