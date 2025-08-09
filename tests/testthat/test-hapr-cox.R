# tests/testthat/test-hapr-cox.R

library(testthat)
library(tidyverse)
library(survival)
library(hapr)

# ---- Helper ----
simulate_mock_dataset <- function(n,
                                  var_v = 1/3,
                                  var_epsilon = 1/3,
                                  var_thetaw = 1/3,
                                  beta_gf = 0.42,
                                  beta_w1 = 0.17,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  w <- data.frame(
    w1 = rnorm(n),
    w2 = factor(sample(c("A", "B", "C"), n, replace = TRUE))
  )

  v <- rnorm(n) * sqrt(var_v)
  epsilon <- rnorm(n) * sqrt(var_epsilon)
  gf <- w$w1 * sqrt(var_thetaw) + v
  gc <- gf + epsilon
  gc_normalized <- scale(gc) |> as.numeric()

  linpred <- beta_gf * gf + beta_w1 * w$w1
  baseline_time <- rexp(n, rate = exp(linpred))
  censoring <- rexp(n, rate = 1)
  time <- pmin(baseline_time, censoring)
  event <- as.integer(baseline_time <= censoring)

  y <- Surv(time, event)

  list(
    y = y,
    gf = gf,
    gc = gc_normalized,
    w = w,
    true_improvement_ratio = 1 / (1 - var_epsilon),
    beta_w1 = beta_w1,
    beta_gf = beta_gf
  )
}

# ---- CI COVERAGE TEST ----
test_that("hapr confidence intervals for beta achieve exact match coverage (cox)", {
  skip_if_not(Sys.getenv("HAPR_RUN_SLOW_TESTS") == "true",
              "Skipping slow coverage test; set HAPR_RUN_SLOW_TESTS=true to run")
  default_params <- list(
    n = 1e4,
    var_v = 1 / 3,
    var_epsilon = 1 / 3,
    var_thetaw = 1 / 3,
    beta_gf = 0.42,
    beta_w1 = 0.17,
    seed = 123
  )

  create_simulated_dataset <- function(params = list()) {
    p <- modifyList(default_params, params)
    set.seed(p$seed)

    n <- p$n
    var_v <- p$var_v
    var_epsilon <- p$var_epsilon
    var_thetaw <- p$var_thetaw
    beta_gf <- p$beta_gf
    beta_w1 <- p$beta_w1

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

    linpred <- beta_gf * gf + beta_w1 * w$w1
    baseline_time <- rexp(n, rate = exp(linpred))
    censoring <- rexp(n, rate = 1)
    time <- pmin(baseline_time, censoring)
    event <- as.integer(baseline_time <= censoring)

    y <- Surv(time, event)

    list(
      y = y,
      gc = gc_normalized,
      w = w,
      true_improvement_ratio = true_improvement_ratio
    )
  }

  sim_data <- create_simulated_dataset()
  fit <- hapr(sim_data$y, sim_data$gc, sim_data$w,
              model_type = "cox",
              improvement_ratio = sim_data$true_improvement_ratio)

  beta_names <- names(fit$coefficients$beta)

  true_beta <- c(gf = default_params$beta_gf,
                 w1 = default_params$beta_w1,
                 w2B = 0,
                 w2C = 0)

  n_sim <- 1000
  covered_matrix <- matrix(NA, nrow = n_sim, ncol = length(beta_names))
  colnames(covered_matrix) <- beta_names

  for (i in 1:n_sim) {
    if (i %% 50 == 0) cat("Simulation", i, "\n")

    sim_data_i <- create_simulated_dataset(params = list(seed = i))
    fit_i <- hapr(
      y = sim_data_i$y,
      gc = sim_data_i$gc,
      w = sim_data_i$w,
      model_type = "cox",
      improvement_ratio = sim_data_i$true_improvement_ratio
    )

    beta_hat <- fit_i$coefficients$beta
    ci_beta <- fit_i$ci_beta

    for (term in beta_names) {
      ci_lower <- ci_beta[term, "Lower"]
      ci_upper <- ci_beta[term, "Upper"]
      covered_matrix[i, term] <- (true_beta[term] >= ci_lower) & (true_beta[term] <= ci_upper)
    }
  }

  coverage_df <- tibble(
    Term = beta_names,
    Coverage = round(colMeans(covered_matrix), 3)
  )

  print(coverage_df)

  expect_true(all(coverage_df$Coverage > 0.90),
              info = paste("Coverage below threshold:\n",
                           paste(coverage_df$Term, coverage_df$Coverage, collapse = "\n")))
})

# ---- hapr_survfit tests ----
test_that("hapr_survfit works with basic input and returns valid structure", {
  sim <- simulate_mock_dataset(n = 200)
  fit <- hapr(sim$y, sim$gc, sim$w, model_type = "cox", improvement_ratio = sim$true_improvement_ratio)

  newdata <- cbind(sim$w, gf = sim$gf)
  survfit_out <- hapr_survfit(fit, newdata = newdata)
  expect_s3_class(survfit_out, "hapr_survfit")
  expect_true("time" %in% names(survfit_out))
  expect_true("surv" %in% names(survfit_out))
  expect_true(all(dim(survfit_out$surv)[2] == nrow(newdata)))
})

test_that("start.time produces valid conditional survival curves", {
  sim <- simulate_mock_dataset(n = 200)
  fit <- hapr(sim$y, sim$gc, sim$w, model_type = "cox", improvement_ratio = sim$true_improvement_ratio)

  newdata <- cbind(sim$w, gf = sim$gf)
  survfit_cond <- hapr_survfit(fit, newdata = newdata[1:50, ], start.time = 2)

  expect_true(min(survfit_cond$time) >= 2)

  # Use closest available time to 2
  start_idx <- which.min(abs(survfit_cond$time - 2))
  expect_true(all(abs(survfit_cond$surv[start_idx, ] - 1) < 1e-6))
})

test_that("aggregation computes average curve and CIs", {
  sim <- simulate_mock_dataset(n = 300)
  fit <- hapr(sim$y, sim$gc, sim$w, model_type = "cox", improvement_ratio = sim$true_improvement_ratio)

  newdata <- cbind(sim$w, gf = sim$gf)
  pred <- predict(fit, newdata = newdata, covariates = "gf_w", type = "risk")$y_hat_gf_w
  threshold <- quantile(pred, 0.9)

  sf <- hapr_survfit(fit, newdata = newdata, aggregate = pred > threshold,
                     conf.int = TRUE, n.boot = 30)

  expect_true("surv_avg" %in% names(sf))
  expect_true("surv_avg_lower" %in% names(sf))
  expect_true("surv_avg_upper" %in% names(sf))
  expect_true(all(sf$surv_avg_upper >= sf$surv_avg))
  expect_true(all(sf$surv_avg >= sf$surv_avg_lower))
})
