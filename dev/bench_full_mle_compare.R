# Benchmark and compare the full MLE fit using original vs. optimized gradient
# for Weibull survival models.

library(hapr)
if (!"testthat" %in% .packages()) {
  library(testthat)
}
if (!"microbenchmark" %in% .packages()) {
  library(microbenchmark)
}

source("tests/testthat/helper-mock_dataset.R")

set.seed(123)

# --- 1. Data Generation ---
n <- 1000 # Using a smaller n for full MLE benchmark to avoid very long run times
p <- 5

beta_g <- 0.6
log_k <- 0.5
censor_rate <- 0.2
var_epsilon_true <- 0.7
var_v_true <- (1 - var_epsilon_true) * 0.4
improvement_ratio_true <- 1 / (1 - var_epsilon_true)

beta_w <- c(0.1, rep(0.05, p))
theta <- c(0.0, rep(0.1, p))

mock_data <- mock_dataset_survival_weibull(
  n = n,
  var_v = var_v_true,
  var_epsilon = var_epsilon_true,
  beta_g = beta_g,
  beta_w = beta_w,
  theta = theta,
  log_k = log_k,
  censor_rate = censor_rate
)

# Initial parameters for optimization
start_beta <- c(0.4, -0.4, 0.05, -0.15, 0.25, -0.05, 0.02) # Adjusted based on p=5
start_delta <- c(log_k = 0)

# --- 2. Test for Equality ---
test_that("Full MLE with optimized and original gradients return identical results", {
  message("Running full MLE with original gradient...")
  original_fit <- hapr_mle_survival(
    event_time = mock_data$event_time,
    event_status = mock_data$event_status,
    gc = mock_data$gc,
    w = mock_data$w,
    improvement_ratio = improvement_ratio_true,
    model_type = "weibull",
    start_beta = start_beta,
    start_delta = start_delta,
    use_analytic_gradient = TRUE,
    use_openmp = FALSE,
    use_optimized_gradient = FALSE, # Use original gradient
    control = list(maxit = 150, reltol = 1e-10, abstol = 1e-10) # Stricter tolerance for equality test
  )

  message("Running full MLE with optimized gradient...")
  optimized_fit <- hapr_mle_survival(
    event_time = mock_data$event_time,
    event_status = mock_data$event_status,
    gc = mock_data$gc,
    w = mock_data$w,
    improvement_ratio = improvement_ratio_true,
    model_type = "weibull",
    start_beta = start_beta,
    start_delta = start_delta,
    use_analytic_gradient = TRUE,
    use_openmp = FALSE,
    use_optimized_gradient = TRUE, # Use optimized gradient
    control = list(maxit = 150, reltol = 1e-10, abstol = 1e-10) # Stricter tolerance for equality test
  )

  # Compare estimated parameters
  expect_equal(original_fit$parameters$beta, optimized_fit$parameters$beta, tolerance = 1e-6)
  expect_equal(original_fit$parameters$delta, optimized_fit$parameters$delta, tolerance = 1e-6)
  expect_equal(original_fit$opt$value, optimized_fit$opt$value, tolerance = 1e-6)
  
  # Compare convergence status
  expect_equal(original_fit$opt$convergence, optimized_fit$opt$convergence)
  expect_equal(original_fit$opt$message, optimized_fit$opt$message)
})

cat("✅ Equality tests passed.\n\n")

# --- 3. Benchmark Performance ---
cat("Running full MLE benchmark...\n")

benchmark_results <- microbenchmark(
  original_gradient = {
    hapr_mle_survival(
      event_time = mock_data$event_time,
      event_status = mock_data$event_status,
      gc = mock_data$gc,
      w = mock_data$w,
      improvement_ratio = improvement_ratio_true,
      model_type = "weibull",
      start_beta = start_beta,
      start_delta = start_delta,
      use_analytic_gradient = TRUE,
      use_openmp = FALSE,
      use_optimized_gradient = FALSE,
      control = list(maxit = 150)
    )
  },
  optimized_gradient = {
    hapr_mle_survival(
      event_time = mock_data$event_time,
      event_status = mock_data$event_status,
      gc = mock_data$gc,
      w = mock_data$w,
      improvement_ratio = improvement_ratio_true,
      model_type = "weibull",
      start_beta = start_beta,
      start_delta = start_delta,
      use_analytic_gradient = TRUE,
      use_openmp = FALSE,
      use_optimized_gradient = TRUE,
      control = list(maxit = 150)
    )
  },
  times = 5L # Fewer repetitions for full MLE due to longer runtime
)

print(benchmark_results)

cat("\nFull MLE benchmark complete. The optimized gradient implementation is ready for review.\n")
