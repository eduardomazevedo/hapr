# Benchmark and compare the original vs. optimized gradient likelihood calculation
# for survival models.

if (!"devtools" %in% .packages()) {
  library(devtools)
}
devtools::load_all()
if (!"Rcpp" %in% .packages()) {
  library(Rcpp)
}
if (!"testthat" %in% .packages()) {
  library(testthat)
}
if (!"microbenchmark" %in% .packages()) {
  library(microbenchmark)
}




# Source the helper file for mock data generation
source("tests/testthat/helper-mock_dataset.R")

# --- 1. Data Generation ---
set.seed(123)
n <- 25000
p <- 5
mock_data <- mock_dataset_survival_weibull(
  n = n,
  var_v = 0.6,
  var_epsilon = 0.1,
  beta_g = 0.5,
  beta_w = c(-0.5, 0.1, -0.2, 0.3, -0.1, 0.05), # Intercept + 5 Wariables
  theta = c(0, 0.5, -0.5, 0.2, -0.2, 0.1), # Intercept + 5 Wariables
  log_k = log(1.5),
  censor_rate = 0.1
)

# --- 2. Prepare Inputs for Likelihood ---
# Run stage 1 to get the conditional expectation of gf
stage1_results <- hapr_first_stage(
  y = mock_data$event_time,
  gc = mock_data$gc,
  w = mock_data$w,
  model_type = "mle"
)

# --- 2b. Calculate avg_linpred (was missing) ---
# This logic is adapted from hapr_mle_survival
improvement_ratio <- stage1_results$stats$max_improvement_ratio / 2
var_epsilon <- 1 - 1 / improvement_ratio
var_v_plus_var_epsilon <- stage1_results$parameters$var_v_plus_var_epsilon
var_v <- var_v_plus_var_epsilon - var_epsilon
if (var_v <= 0) {
  stop("Derived var_v is not positive, check improvement_ratio")
}
posterior <- abc(var_epsilon, var_v)
theta_hat <- stage1_results$parameters$theta
X_w_int <- cbind(1, stage1_results$preprocessed$w)
colnames(X_w_int)[1] <- "(Intercept)"
w_theta <- c(X_w_int %*% theta_hat)
avg_linpred <- posterior$a * stage1_results$preprocessed$gc + posterior$b * w_theta

# Create design matrices
X <- cbind(1, mock_data$w)
event_indices <- which(mock_data$event_status == 1)
censor_indices <- which(mock_data$event_status == 0)

# Event data
event_time_event <- mock_data$event_time[event_indices]
avg_linpred_event <- avg_linpred[event_indices]
X_w_event <- X[event_indices, ]

# Censored data
event_time_censor <- mock_data$event_time[censor_indices]
avg_linpred_censor <- avg_linpred[censor_indices]
X_w_censor <- X[censor_indices, ]

# Other parameters
post_c <- posterior$c
model_type <- 1 # Weibull

# Starting parameters for the likelihood function (beta_g, beta_w, log_k)
params <- c(0.4, -0.4, 0.05, -0.15, 0.25, -0.05, 0.02, log(1.6))

# --- 3. Test for Equality ---
test_that("Optimized and original functions return identical results", {
  # Run original function
  original_res <- hapr_mle_survival_nll_split_grad_cpp(
    params = params,
    event_time = event_time_event,
    avg_linpred_event = avg_linpred_event,
    X_w_event = X_w_event,
    censor_time = event_time_censor,
    avg_linpred_censor = avg_linpred_censor,
    X_w_censor = X_w_censor,
    post_c = post_c,
    model_type = model_type,
    use_openmp = FALSE # Use single thread for fair comparison
  )

  # Run optimized function
  optimized_res <- hapr_mle_survival_nll_split_grad_cpp_optimized(
    params = params,
    event_time = event_time_event,
    avg_linpred_event = avg_linpred_event,
    X_w_event = X_w_event,
    censor_time = event_time_censor,
    avg_linpred_censor = avg_linpred_censor,
    X_w_censor = X_w_censor,
    post_c = post_c,
    model_type = model_type,
    use_openmp = FALSE
  )

  # Compare value
  expect_equal(original_res$value, optimized_res$value, tolerance = 1e-9)

  # Compare gradient
  expect_equal(original_res$gradient, optimized_res$gradient, tolerance = 1e-9)
})

cat("✅ Equality tests passed.\n\n")


# --- 4. Benchmark Performance ---
cat("Running benchmark...\n")
benchmark_results <- microbenchmark(
  original = {
    hapr_mle_survival_nll_split_grad_cpp(
      params = params,
      event_time = event_time_event,
      avg_linpred_event = avg_linpred_event,
      X_w_event = X_w_event,
      censor_time = event_time_censor,
      avg_linpred_censor = avg_linpred_censor,
      X_w_censor = X_w_censor,
      post_c = post_c,
      model_type = model_type,
      use_openmp = FALSE
    )
  },
  optimized = {
    hapr_mle_survival_nll_split_grad_cpp_optimized(
      params = params,
      event_time = event_time_event,
      avg_linpred_event = avg_linpred_event,
      X_w_event = X_w_event,
      censor_time = event_time_censor,
      avg_linpred_censor = avg_linpred_censor,
      X_w_censor = X_w_censor,
      post_c = post_c,
      model_type = model_type,
      use_openmp = FALSE
    )
  },
  times = 20L # Number of repetitions
)

print(benchmark_results)

cat("\nBenchmark complete. The new optimized implementation is ready for review.\n")
