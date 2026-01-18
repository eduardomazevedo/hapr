#' Script to test hapr_mle_survival with Weibull survival data
#'
#' Creates a dataset with 3 covariates and n=1000, then runs Weibull MLE.

devtools::load_all()

source("tests/testthat/helper-mock_dataset.R")

set.seed(123)

n <- 1000
var_epsilon <- 0.7
var_v <- (1 - var_epsilon) * 0.5

beta_g <- 0.6
beta_w <- c(0.1, -0.2, 0.15, 0.05)
theta <- c(0.0, 0.1, -0.25, 0.2)
log_k <- 0

improvement_ratio <- 1 / (1 - var_epsilon)

cat("Creating mock dataset (weibull)...\n")
data <- mock_dataset_survival_weibull(
  n = n,
  var_v = var_v,
  var_epsilon = var_epsilon,
  beta_g = beta_g,
  beta_w = beta_w,
  theta = theta,
  log_k = log_k,
  censor_rate = 0.2
)

cat(sprintf("Dataset created: n=%d, %d covariates\n", n, ncol(data$w)))
cat("\n")

start_beta <- rep(0, ncol(data$w) + 2)
start_delta <- c(log_k = 0)

cat("Running hapr_mle_survival (weibull)...\n")
mle_fit <- hapr_mle_survival(
  event_time = data$event_time,
  event_status = data$event_status,
  gc = data$gc,
  w = data$w,
  improvement_ratio = improvement_ratio,
  model_type = "weibull",
  start_beta = start_beta,
  start_delta = start_delta,
  control = list(maxit = 150)
)

print(mle_fit)
