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
n <- 100
var_v <- 0.3
var_epsilon <- 0.5
beta_g <- 0.42  # Effect of future PRS (gf) on outcome
beta_w <- c(0.1, 0.17, 0.27, -0.27)  # Intercept + 3 covariate effects
theta <- c(0.0, 0.1, -0.2, 0.3)  # Intercept + 3 covariate effects for gc ~ w
var_y <- 1.0  # Error variance for outcome

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

# Print results
cat("\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("HAPR First Stage Results\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

cat("PARAMETERS:\n")
cat("Theta (gc ~ w):\n")
print(first_stage_fit$parameters$theta)
cat("\n")

cat("Var(v + epsilon):\n")
print(first_stage_fit$parameters$var_v_plus_var_epsilon)
cat("\n")

cat("Gamma (y ~ gc + w):\n")
print(first_stage_fit$parameters$gamma)
cat("\n")

cat("STATISTICS:\n")
cat("Max improvement ratio:", first_stage_fit$stats$max_improvement_ratio, "\n")
cat("Var(w*theta):", first_stage_fit$stats$var_wtheta, "\n")
cat("\n")

cat("REGRESSION R-SQUARED:\n")
cat("gc ~ w R²:", first_stage_fit$regressions$gc_on_w$r2, "\n")
cat("y ~ w R²:", first_stage_fit$regressions$y_on_w$r2, "\n")
cat("y ~ gc R²:", first_stage_fit$regressions$y_on_gc$r2, "\n")
cat("y ~ gc + w R²:", first_stage_fit$regressions$y_on_gc_w$r2, "\n")
cat("\n")

cat(paste(rep("=", 72), collapse = ""), "\n")

print(first_stage_fit)
