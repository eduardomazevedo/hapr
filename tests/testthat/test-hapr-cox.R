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
    gc = gc_normalized,
    w = w,
    true_improvement_ratio = 1 / (1 - var_epsilon),
    beta_w1 = beta_w1,
    beta_gf = beta_gf
  )
}

# ---- CI COVERAGE TEST ----
test_that("hapr confidence intervals for beta achieve exact match coverage (cox)", {
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
  true_beta <- fit$coefficients$beta

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

  expect_true(all(coverage_df$Coverage > 0.70),
              info = paste("Coverage below threshold:\n",
                           paste(coverage_df$Term, coverage_df$Coverage, collapse = "\n")))
})
