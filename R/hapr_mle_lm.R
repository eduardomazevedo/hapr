#' HAPR maximum likelihood estimation
#'
#' @description
#' Fits the HAPR linear model by maximum likelihood using Gauss-Hermite quadrature.
#'
#' @param y Outcome vector
#' @param gc Current PRS vector (will be normalized)
#' @param w Covariate matrix (no intercept column)
#' @param improvement_ratio R-squared improvement ratio (required)
#' @param start_beta Named numeric vector for beta parameters (gf, (Intercept), w1, ...)
#' @param start_delta Named numeric vector for additional parameters (optional, log_sigma)
#' @param control List passed to stats::optim
#'
#' @return A hapr_mle_fit object with MLE estimates and diagnostics
#' @export
hapr_mle <- function(
    y,
    gc,
    w,
    improvement_ratio,
    start_beta,
    start_delta = NULL,
    control = list()) {
  if (missing(improvement_ratio) || is.null(improvement_ratio)) {
    stop("improvement_ratio must be specified.")
  }
  if (is.null(start_beta) || !is.numeric(start_beta)) {
    stop("start_beta must be a numeric vector")
  }
  if (any(is.na(y))) {
    stop("y must not contain missing values")
  }

  first_stage <- hapr_first_stage(y = y, gc = gc, w = w, model_type = "mle")

  gc <- first_stage$preprocessed$gc
  w <- first_stage$preprocessed$w

  add_intercept <- function(X) {
    X_with_int <- cbind(1, X)
    colnames(X_with_int)[1] <- "(Intercept)"
    X_with_int
  }

  theta_hat <- first_stage$parameters$theta
  var_v_plus_var_epsilon <- first_stage$parameters$var_v_plus_var_epsilon
  max_improvement_ratio <- first_stage$stats$max_improvement_ratio
  if (improvement_ratio >= max_improvement_ratio) {
    stop(sprintf("Improvement ratio must be less than %s.", max_improvement_ratio))
  }

  var_epsilon <- 1 - 1 / improvement_ratio
  var_v <- var_v_plus_var_epsilon - var_epsilon
  if (var_v <= 0) {
    stop("Derived var_v must be positive; check improvement_ratio.")
  }
  posterior <- abc(var_epsilon, var_v)

  # Prepare beta parameter names.
  nb <- ncol(w) + 2
  beta_names <- c("gf", "(Intercept)", colnames(w))
  if (length(start_beta) != nb) {
    stop(sprintf("start_beta must have length %d.", nb))
  }
  names(start_beta) <- beta_names

  if (!is.null(start_delta)) {
    if (!is.numeric(start_delta)) {
      stop("start_delta must be a numeric vector")
    }
    if (length(start_delta) != 1) {
      stop("start_delta must contain exactly one element (log_sigma).")
    }
    names(start_delta) <- "log_sigma"
  } else {
    start_delta <- c(log_sigma = 0)
  }

  start_params <- c(start_beta, start_delta)

  X_w <- add_intercept(w)
  w_theta <- c(X_w %*% theta_hat)
  neg_loglik <- make_hapr_mle_likelihood_lm(
    y = y,
    gc = gc,
    w_theta = w_theta,
    X_w = X_w,
    posterior = posterior
  )
  
  opt <- stats::optim(
    par = start_params,
    fn = neg_loglik,
    method = "BFGS",
    control = control,
    hessian = TRUE
  )

  mle_params <- opt$par
  beta_hat <- mle_params[seq_len(nb)]
  names(beta_hat) <- beta_names
  if (length(mle_params) > nb) {
    delta_hat <- mle_params[(nb + 1):length(mle_params)]
    names(delta_hat) <- names(start_delta)
  } else {
    delta_hat <- numeric(0)
  }

  vcov_all <- NULL
  standard_errors <- NULL
  ci_beta <- NULL
  if (is.matrix(opt$hessian)) {
    vcov_try <- try(solve(opt$hessian), silent = TRUE)
    if (!inherits(vcov_try, "try-error")) {
      dimnames(vcov_try) <- list(names(mle_params), names(mle_params))
      vcov_all <- vcov_try
      standard_errors <- sqrt(diag(vcov_all))[seq_len(nb)]
      names(standard_errors) <- beta_names
      z <- stats::qnorm(0.975)
      ci_beta <- data.frame(
        Estimate = beta_hat,
        Std.Error = standard_errors,
        Lower = beta_hat - z * standard_errors,
        Upper = beta_hat + z * standard_errors,
        row.names = beta_names,
        check.names = FALSE
      )
    }
  }

  result <- list(
    model_type = "mle",
    regressions = list(gc_on_w = first_stage$regressions$gc_on_w),
    parameters = list(
      beta = beta_hat,
      delta = delta_hat,
      theta = theta_hat,
      var_v_plus_var_epsilon = var_v_plus_var_epsilon
    ),
    vcov_parameters = list(
      all = vcov_all
    ),
    standard_errors = standard_errors,
    ci_beta = ci_beta,
    stats = list(
      var_v = var_v,
      var_epsilon = var_epsilon,
      posterior = posterior,
      improvement_ratio = improvement_ratio,
      max_improvement_ratio = max_improvement_ratio
    ),
    opt = opt
  )
  class(result) <- "hapr_mle_fit"
  result
}
