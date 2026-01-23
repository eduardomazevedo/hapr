#' Benchmark Weibull survival MLE across n and p

suppressMessages({
  devtools::load_all()
})

source("tests/testthat/helper-mock_dataset.R")

set.seed(123)

ns <- c(10000, 50000, 100000, 250000, 500000)
ps <- c(2, 5, 10, 20, 30)

candidates <- expand.grid(n = ns, p = ps, KEEP.OUT.ATTRS = FALSE)
# Keep at most 10 combos by spreading across n and p
pick_idx <- unique(round(seq(1, nrow(candidates), length.out = 10)))
combos <- candidates[pick_idx, , drop = FALSE]

var_epsilon <- 0.7
var_v <- (1 - var_epsilon) * 0.4
improvement_ratio <- 1 / (1 - var_epsilon)

beta_g <- 0.6
log_k <- 0.5
censor_rate <- 0.2

benchmark_once <- function(n, p, use_openmp) {
  beta_w <- c(0.1, rep(0.05, p))
  theta <- c(0.0, rep(0.1, p))

  data <- mock_dataset_survival_weibull(
    n = n,
    var_v = var_v,
    var_epsilon = var_epsilon,
    beta_g = beta_g,
    beta_w = beta_w,
    theta = theta,
    log_k = log_k,
    censor_rate = censor_rate
  )

  start_beta <- rep(0, ncol(data$w) + 2)

  elapsed <- system.time({
    hapr_mle_survival(
      event_time = data$event_time,
      event_status = data$event_status,
      gc = data$gc,
      w = data$w,
      improvement_ratio = improvement_ratio,
      model_type = "weibull",
      start_beta = start_beta,
      start_delta = c(log_k = 0),
      use_openmp = use_openmp,
      control = list(maxit = 150)
    )
  })[["elapsed"]]

  elapsed
}

results <- list()

for (i in seq_len(nrow(combos))) {
  n <- combos$n[i]
  p <- combos$p[i]

  elapsed_time <- benchmark_once(n, p, use_openmp = TRUE)

  results[[i]] <- data.frame(
    n = n,
    p = p,
    elapsed_s = elapsed_time,
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, results)
summary_table <- summary_table[order(summary_table$n, summary_table$p), , drop = FALSE]

out_dir <- "dev/output"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_file <- file.path(out_dir, "bench_weibull_mle.csv")
write.csv(summary_table, out_file, row.names = FALSE)

print(summary_table)
cat(sprintf("\nWrote %s\n", out_file))
