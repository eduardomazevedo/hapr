#' HAPR maximum likelihood estimation for survival models
#'
#' @description
#' Fits the HAPR parametric survival models by maximum likelihood using
#' Gauss-Hermite quadrature. Supported models are exponential and Weibull.
#'
#' @param event_time Event or censoring time (numeric vector)
#' @param event_status Event indicator (0/1 or logical)
#' @param gc Current PRS vector (will be normalized)
#' @param w Covariate matrix (no intercept column)
#' @param improvement_ratio R-squared improvement ratio (required)
#' @param model_type One of "exponential" or "weibull"
#' @param start_beta Named numeric vector for beta parameters (gf, (Intercept), w1, ...)
#' @param start_delta Named numeric vector for additional parameters (optional);
#'   for Weibull this should contain log_k.
#' @param control List passed to stats::optim
#'
#' @return A hapr_mle_fit object with MLE estimates and diagnostics
#' @export
hapr_mle_survival <- function(
    event_time,
    event_status,
    gc,
    w,
    improvement_ratio,
    model_type = "exponential",
    start_beta,
    start_delta = NULL,
    control = list()) {
  model_type <- match.arg(model_type, c("exponential", "weibull"))
  if (missing(improvement_ratio) || is.null(improvement_ratio)) {
    stop("improvement_ratio must be specified.")
  }
  if (is.null(start_beta) || !is.numeric(start_beta)) {
    stop("start_beta must be a numeric vector")
  }
  if (!is.numeric(event_time)) {
    stop("event_time must be numeric")
  }
  if (any(is.na(event_time))) {
    stop("event_time must not contain missing values")
  }
  if (any(event_time < 0)) {
    stop("event_time must be non-negative")
  }
  if (is.logical(event_status)) {
    event_status <- as.numeric(event_status)
  }
  if (!is.numeric(event_status)) {
    stop("event_status must be numeric or logical")
  }
  if (any(is.na(event_status))) {
    stop("event_status must not contain missing values")
  }
  if (any(!event_status %in% c(0, 1))) {
    stop("event_status must be 0/1 or logical")
  }
  if (length(event_time) != length(event_status)) {
    stop("event_time and event_status must have the same length")
  }

  first_stage <- hapr_first_stage(
    y = event_time,
    gc = gc,
    w = w,
    model_type = "mle"
  )

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
    if (model_type == "exponential" && length(start_delta) != 0) {
      stop("start_delta must be empty for exponential survival.")
    }
    if (model_type == "weibull" && length(start_delta) != 1) {
      stop("start_delta must contain exactly one element (log_k).")
    }
  } else if (model_type == "weibull") {
    start_delta <- c(log_k = 0)
  } else {
    start_delta <- numeric(0)
  }
  if (model_type == "weibull") {
    names(start_delta) <- "log_k"
  }

  start_params <- c(start_beta, start_delta)

  X_w <- add_intercept(w)
  w_theta <- c(X_w %*% theta_hat)
  neg_loglik <- make_hapr_mle_likelihood_survival(
    event_time = event_time,
    event_status = event_status,
    gc = gc,
    w_theta = w_theta,
    X_w = X_w,
    posterior = posterior,
    model_type = model_type
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
    model_type = sprintf("mle_survival_%s", model_type),
    regressions = list(gc_on_w = first_stage$regressions$gc_on_w),
    parameters = list(
      beta = beta_hat,
      delta = delta_hat,
      theta = theta_hat,
      var_v_plus_var_epsilon = var_v_plus_var_epsilon
    ),
    vcov_parameters = list(
      all = vcov_all,
      order = list(
        beta = seq_len(nb),
        delta = if (length(mle_params) > nb) (nb + 1):length(mle_params) else integer(0)
      )
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

make_hapr_mle_likelihood_survival <- function(
    event_time,
    event_status,
    gc,
    w_theta,
    X_w,
    posterior,
    model_type) {
  avg_linpred <- posterior$a * gc + posterior$b * w_theta
  event_idx <- which(event_status == 1)
  censor_idx <- which(event_status == 0)

  event_time_event <- event_time[event_idx]
  avg_event <- avg_linpred[event_idx]
  X_w_event <- X_w[event_idx, , drop = FALSE]

  censor_time <- event_time[censor_idx]
  avg_censor <- avg_linpred[censor_idx]
  X_w_censor <- X_w[censor_idx, , drop = FALSE]

  if (length(event_idx) == 0) {
    X_w_event <- X_w[0, , drop = FALSE]
  }
  if (length(censor_idx) == 0) {
    X_w_censor <- X_w[0, , drop = FALSE]
  }
  model_id <- if (model_type == "exponential") 0L else 1L
  function(params) {
    hapr_mle_survival_nll_split_cpp(
      params = params,
      event_time = event_time_event,
      avg_linpred_event = avg_event,
      X_w_event = X_w_event,
      censor_time = censor_time,
      avg_linpred_censor = avg_censor,
      X_w_censor = X_w_censor,
      post_c = posterior$c,
      model_type = model_id
    )
  }
}
