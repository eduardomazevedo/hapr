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
    start_beta,
    start_delta = NULL,
    use_openmp = TRUE,
    control = list()) {
  model_type <- match.arg(model_type, c("exponential", "weibull"))
  if (missing(improvement_ratio) || is.null(improvement_ratio)) {
    stop("improvement_ratio must be specified.")
  }
  if (!is.logical(use_openmp) || length(use_openmp) != 1) {
    stop("use_openmp must be a single logical value.")
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

  # Convert start_beta -> start_gamma for gamma-parameterized likelihood
  start_gamma <- start_beta
  start_gamma["gf"] <- start_beta["gf"] * posterior$a
  if (nb > 1) {
    start_gamma[-1] <- start_beta[-1] + start_beta["gf"] * posterior$b * theta_hat
  }
  names(start_gamma)[1] <- "gc"

  hapr_mle_survival_gamma(
    event_time = event_time,
    event_status = event_status,
    gc = gc,
    w = w,
    improvement_ratio = improvement_ratio,
    model_type = model_type,
    start_gamma = start_gamma,
    start_delta = start_delta,
    use_openmp = use_openmp,
    control = control
  )
}
