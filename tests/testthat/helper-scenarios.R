stage2_scenarios <- function(run_slow,
                             include_large_n = TRUE,
                             var_v_factors = c(0.5, 0.9, 0.99)) {
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
        for (var_v_factor in var_v_factors) {
          idx_var_v_factor <- which(var_v_factors == var_v_factor)
          scenario_name <- sprintf(
            "%s_n%d_ve%.1f_vv%.2f",
            model_type, as.integer(n), var_epsilon, var_v_factor
          )
          seed <- 123 +
            idx_var_epsilon * 1000 +
            idx_n * 100 +
            idx_model_type * 10 +
            idx_var_v_factor

          var_v <- (1 - var_epsilon) * var_v_factor
          improvement_ratio <- 1 / (1 - var_epsilon)

          scenarios[[scenario_name]] <- list(
            name = scenario_name,
            var_epsilon = var_epsilon,
            n = n,
            model_type = model_type,
            idx_var_epsilon = idx_var_epsilon,
            idx_n = idx_n,
            idx_model_type = idx_model_type,
            idx_var_v_factor = idx_var_v_factor,
            var_v_factor = var_v_factor,
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

#' Build survival test scenarios across model types and var_v factors.
#'
#' @param run_slow Logical. Whether to include large n values.
#' @param log_k_values Numeric vector of log_k values for Weibull scenarios.
#' @param include_large_n Logical. Include n = 1e5 when run_slow is TRUE.
#' @param var_v_factors Numeric vector. Factors for var_v as a fraction of (1 - var_epsilon).
#'
#' @return Named list of scenario definitions.
survival_scenarios <- function(run_slow,
                               log_k_values,
                               include_large_n = TRUE,
                               var_v_factors = c(0.5, 0.9, 0.99)) {
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
          for (var_v_factor in var_v_factors) {
            idx_var_v_factor <- which(var_v_factors == var_v_factor)
            scenario_name <- sprintf("%s_n%d_ve%.1f_vv%.2f_lk%.1f",
                                     substr(model_type, 1, 3),
                                     as.integer(n),
                                     var_epsilon,
                                     var_v_factor,
                                     log_k)

            seed <- 123 +
              idx_var_epsilon * 1000 +
              idx_n * 100 +
              idx_model_type * 10 +
              idx_var_v_factor * 2 +
              round(log_k * 10)

            var_v <- (1 - var_epsilon) * var_v_factor
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
              idx_var_v_factor = idx_var_v_factor,
              var_v_factor = var_v_factor,
              var_v = var_v,
              improvement_ratio = improvement_ratio,
              seed = seed
            )
          }
        }
      }
    }
  }

  scenarios
}
