#' Script to test hapr_first_stage with mock data
#' 
#' Creates a dataset with 3 covariates and n=100, then runs hapr stage 1

# Load package (run from project root)
devtools::load_all()

# Source mock dataset functions
source("tests/testthat/helper-mock_dataset.R")

# Set seed for reproducibility
set.seed(123)

# Parameters for mock dataset
n <- 1000
var_epsilon <- 0.9
var_v <- (1 - var_epsilon) * 0.5

beta_g <- 1.42  # Effect of future PRS (gf) on outcome
beta_w <- c(0.1, 0.17, 0.27, -0.27)  # Intercept + 3 covariate effects
theta <- c(0.0, 0.1, -0.2, 0.3)  # Intercept + 3 covariate effects for gc ~ w
var_y <- 1.0  # Error variance for outcome

improvement_ratio <- 1 / (1 - var_epsilon)

# Create dataset
cat("Creating mock dataset...\n")
data <- mock_dataset_lm(
  n = n,
  var_v = var_v,
  var_epsilon = var_epsilon,
  beta_g = beta_g,
  beta_w = beta_w,
  theta = theta,
  var_y = var_y
)

# Extract components
w <- data$w
gc <- data$gc
y <- data$y

cat(sprintf("Dataset created: n=%d, %d covariates\n", n, ncol(w)))
cat("\n")

# Run hapr first stage
cat("Running hapr_first_stage...\n")
first_stage_fit <- hapr_first_stage(
  y = y,
  gc = gc,
  w = w,
  model_type = "lm"
)

# Second stage
cat("Running hapr_second_stage...\n")
second_stage_fit <- hapr_second_stage(
  first_stage = first_stage_fit,
  improvement_ratio = improvement_ratio
)

print(second_stage_fit)
