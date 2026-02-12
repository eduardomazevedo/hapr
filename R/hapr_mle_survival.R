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
#' @param start_beta Optional named numeric vector for beta parameters
#'   (gf, (Intercept), w1, ...). If omitted, a data-driven initialization is used.
#' @param start_delta Named numeric vector for additional parameters (optional);
#'   for Weibull this should contain log_k.
#' @param start_gamma_method Method for initializing gamma parameters.
#'   One of "auto", "beta_transform", or "log_time_lm".
#' @param use_openmp Logical; if TRUE uses OpenMP in the analytic gradient.
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
    start_beta = NULL,
    start_delta = NULL,
    start_gamma_method = c("auto", "beta_transform", "log_time_lm"),
    use_openmp = TRUE,
    control = list()) {
  model_type <- match.arg(model_type, c("exponential", "weibull"))
  start_gamma_method <- match.arg(start_gamma_method)
  if (missing(improvement_ratio) || is.null(improvement_ratio)) {
    stop("improvement_ratio must be specified.")
  }
  if (!is.logical(use_openmp) || length(use_openmp) != 1) {
    stop("use_openmp must be a single logical value.")
  }
  if (!is.null(start_beta) && !is.numeric(start_beta)) {
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
  if (is.null(start_beta)) {
    start_beta <- rep(0, nb)
  } else if (length(start_beta) != nb) {
    stop(sprintf("start_beta must have length %d when provided.", nb))
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

  get_start_gamma_from_beta <- function() {
    start_gamma <- start_beta
    start_gamma["gf"] <- start_beta["gf"] * posterior$a
    if (nb > 1) {
      start_gamma[-1] <- start_beta[-1] + start_beta["gf"] * posterior$b * theta_hat
    }
    names(start_gamma)[1] <- "gc"
    start_gamma
  }

  get_start_gamma_from_log_time <- function() {
    signed_log_time <- log(pmax(event_time, .Machine$double.xmin))
    if (model_type == "exponential") {
      signed_log_time <- -signed_log_time
    }
    x_gc_w <- cbind(gc = gc, `(Intercept)` = 1, w)
    fit <- stats::lm.fit(x = x_gc_w, y = signed_log_time)
    coef <- fit$coefficients
    if (length(coef) != nb) {
      return(NULL)
    }
    coef[is.na(coef)] <- 0
    if (any(!is.finite(coef))) {
      return(NULL)
    }
    names(coef) <- c("gc", "(Intercept)", colnames(w))
    coef
  }

  start_gamma <- switch(
    start_gamma_method,
    beta_transform = get_start_gamma_from_beta(),
    log_time_lm = {
      log_time_start <- get_start_gamma_from_log_time()
      if (is.null(log_time_start)) get_start_gamma_from_beta() else log_time_start
    },
    auto = {
      log_time_start <- get_start_gamma_from_log_time()
      if (is.null(log_time_start)) get_start_gamma_from_beta() else log_time_start
    }
  )

  start_params <- c(start_gamma, start_delta)

  X_w <- add_intercept(w)
  neg_loglik <- make_hapr_mle_likelihood_survival_grad(
    event_time = event_time,
    event_status = event_status,
    gc = gc,
    X_w = X_w,
    posterior = posterior,
    model_type = model_type,
    use_openmp = use_openmp
  )

  opt <- stats::optim(
    par = start_params,
    fn = neg_loglik$fn,
    gr = neg_loglik$gr,
    method = "BFGS",
    control = control,
    hessian = TRUE
  )

  mle_params <- opt$par
  gamma_hat <- mle_params[seq_len(nb)]
  gamma_names <- c("gc", "(Intercept)", colnames(w))
  names(gamma_hat) <- gamma_names
  if (length(mle_params) > nb) {
    delta_hat <- as.list(mle_params[(nb + 1):length(mle_params)])
    names(delta_hat) <- names(start_delta)
  } else {
    delta_hat <- list()
  }

  beta_hat <- gamma_hat
  beta_hat["gc"] <- gamma_hat["gc"] / posterior$a
  if (nb > 1) {
    beta_hat[-1] <- gamma_hat[-1] - beta_hat["gc"] * posterior$b * theta_hat
  }
  names(beta_hat)[1] <- "gf"

  vcov_all <- NULL
  vcov_gamma <- NULL
  standard_errors_gamma <- NULL
  vcov_beta <- NULL
  standard_errors_beta <- NULL
  ci_beta <- NULL

  if (is.matrix(opt$hessian)) {
    vcov_try <- try(solve(opt$hessian), silent = TRUE)
    if (!inherits(vcov_try, "try-error")) {
      dimnames(vcov_try) <- list(names(mle_params), names(mle_params))
      vcov_all <- vcov_try

      vcov_gamma <- vcov_all[seq_len(nb), seq_len(nb), drop = FALSE]
      standard_errors_gamma <- sqrt(diag(vcov_gamma))
      names(standard_errors_gamma) <- gamma_names

      joint_vcov_theta_var_total <- first_stage$vcov_parameters$joint_theta_var_v_plus_var_epsilon
      if (is.null(joint_vcov_theta_var_total)) {
        stop("Missing joint stage-1 covariance for (theta, var_v_plus_var_epsilon).")
      }

      ng <- length(gamma_hat)
      nt <- length(theta_hat)

      vcov_full <- matrix(0, ng + nt + 1, ng + nt + 1)
      vcov_full[1:ng, 1:ng] <- vcov_gamma
      idx_stage1 <- (ng + 1):(ng + nt + 1)
      vcov_full[idx_stage1, idx_stage1] <- joint_vcov_theta_var_total

      J <- calculate_analytical_jacobian(
        model_type = "lm",
        gamma = gamma_hat,
        theta = theta_hat,
        var_total = var_v_plus_var_epsilon,
        posterior = posterior,
        beta = beta_hat,
        derived_vars = list(var_epsilon = var_epsilon)
      )

      vcov_beta <- J %*% vcov_full %*% t(J)
      standard_errors_beta <- sqrt(diag(vcov_beta))
      names(standard_errors_beta) <- names(beta_hat)
      z <- stats::qnorm(0.975)
      ci_beta <- data.frame(
        Estimate = beta_hat,
        Std.Error = standard_errors_beta,
        Lower = beta_hat - z * standard_errors_beta,
        Upper = beta_hat + z * standard_errors_beta,
        row.names = names(beta_hat),
        check.names = FALSE
      )
    }
  }

  result <- list(
    model_type = sprintf("mle_survival_%s", model_type),
    regressions = list(gc_on_w = first_stage$regressions$gc_on_w),
    parameters = list(
      gamma = gamma_hat,
      beta = beta_hat,
      delta = delta_hat,
      theta = theta_hat,
      var_v_plus_var_epsilon = var_v_plus_var_epsilon
    ),
    vcov_parameters = list(
      all = vcov_all,
      order = list(
        gamma = seq_len(nb),
        delta = if (length(mle_params) > nb) (nb + 1):length(mle_params) else integer(0)
      ),
      gamma = vcov_gamma,
      beta = vcov_beta
    ),
    standard_errors = standard_errors_beta,
    standard_errors_gamma = standard_errors_gamma,
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
