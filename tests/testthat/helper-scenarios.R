stage2_scenarios <- function(run_slow, include_large_n = TRUE) {
  var_epsilon_values <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  if (include_large_n && run_slow) {
    n_values <- c(1e3, 1e4, 1e5)
  } else {
    n_values <- c(1e3, 1e4)
  }
  model_types <- c("lm", "probit")

  scenarios <- list()
  for (var_epsilon in var_epsilon_values) {
    idx_var_epsilon <- which(var_epsilon_values == var_epsilon)
    for (n in n_values) {
      idx_n <- which(n_values == n)
      for (model_type in model_types) {
        idx_model_type <- which(model_types == model_type)
        scenario_name <- sprintf("%s_n%d_ve%.1f", model_type, as.integer(n), var_epsilon)
        seed <- 123 +
          idx_var_epsilon * 100 +
          idx_n * 10 +
          idx_model_type

        var_v <- (1 - var_epsilon) * 0.5
        improvement_ratio <- 1 / (1 - var_epsilon)

        scenarios[[scenario_name]] <- list(
          name = scenario_name,
          var_epsilon = var_epsilon,
          n = n,
          model_type = model_type,
          idx_var_epsilon = idx_var_epsilon,
          idx_n = idx_n,
          idx_model_type = idx_model_type,
          var_v = var_v,
          improvement_ratio = improvement_ratio,
          seed = seed
        )
      }
    }
  }

  scenarios
}

survival_scenarios <- function(run_slow, log_k_values, include_large_n = TRUE) {
  var_epsilon_values <- c(0.5, 0.6, 0.7, 0.8, 0.9)
  if (include_large_n && run_slow) {
    n_values <- c(1e3, 1e4, 1e5)
  } else {
    n_values <- c(1e3, 1e4)
  }
  model_types <- c("exponential", "weibull")

  scenarios <- list()
  for (var_epsilon in var_epsilon_values) {
    idx_var_epsilon <- which(var_epsilon_values == var_epsilon)
    for (n in n_values) {
      idx_n <- which(n_values == n)
      for (model_type in model_types) {
        idx_model_type <- which(model_types == model_type)
        model_log_k_values <- if (model_type == "weibull") log_k_values else 0
        for (log_k in model_log_k_values) {
          scenario_name <- sprintf("%s_n%d_ve%.1f_lk%.1f",
                                   substr(model_type, 1, 3),
                                   as.integer(n),
                                   var_epsilon,
                                   log_k)

          seed <- 123 +
            idx_var_epsilon * 100 +
            idx_n * 10 +
            idx_model_type * 1000 +
            round(log_k * 10)

          var_v <- (1 - var_epsilon) * 0.4
          improvement_ratio <- 1 / (1 - var_epsilon)

          scenarios[[scenario_name]] <- list(
            name = scenario_name,
            var_epsilon = var_epsilon,
          n = n,
          model_type = model_type,
          log_k = log_k,
          idx_var_epsilon = idx_var_epsilon,
          idx_n = idx_n,
          idx_model_type = idx_model_type,
          var_v = var_v,
          improvement_ratio = improvement_ratio,
          seed = seed
          )
        }
      }
    }
  }

  scenarios
}
