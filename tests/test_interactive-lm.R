# Create a small interactive test for hapr_lm
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


# Call the hapr_lm function
# fit <- hapr_lm(y, gc_normalized, w, true_improvement_ratio)

first_stage <- hapr_lm_first_stage(y, gc_normalized, w)
second_stage <- hapr_lm_second_stage(first_stage, improvement_ratio = true_improvement_ratio)

# Print the results
# print(fit)
# print(first_stage)
# print(second_stage)
