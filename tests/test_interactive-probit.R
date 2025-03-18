# Create a small interactive test for hapr_probit
set.seed(123)  # For reproducibility
devtools::load_all()

# Create fake data
n <- 1e4

var_v <- 1/3
var_epsilon <- 1/3
var_thetaw <- 1/3

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

y_binary <- as.numeric(y > mean(y)) |> as.factor()

# Call the hapr_probit_first_stage function
first_stage <- hapr_probit_first_stage(y_binary, gc_normalized, w)

# Call the hapr_probit_second_stage function
second_stage <- hapr_probit_second_stage(first_stage, improvement_ratio = 1.5)

# full_fit <- hapr_probit(y_binary, gc_normalized, w, improvement_ratio = 1.5)

# Print the results
print(first_stage)
print(second_stage)
# print(full_fit)
