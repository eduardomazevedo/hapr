#' HAPR maximum likelihood estimation
#'
#' @description
#' Fits the HAPR model by maximum likelihood using a user-supplied log-likelihood
#' function for the outcome conditional on the latent index.
#'
#' @param y Outcome vector
#' @param gc Current PRS vector (will be normalized)
#' @param w Covariate matrix (no intercept column)
#' @param improvement_ratio R-squared improvement ratio (required)
#' @param loglik_fn Function with signature loglik_fn(y, linpred, delta)
#' @param start_beta Named numeric vector for beta parameters (gf, (Intercept), w1, ...)
#' @param start_delta Named numeric vector for additional parameters (optional)
#' @param control List passed to stats::optim
#'
#' @return A hapr_mle_fit object with MLE estimates and diagnostics
#' @export
hapr_mle <- function(
    y,
    gc,
    w,
    improvement_ratio,
    loglik_fn,
    start_beta,
    start_delta = NULL,
    control = list()) {
  if (missing(improvement_ratio) || is.null(improvement_ratio)) {
    stop("improvement_ratio must be specified.")
  }
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function")
  }
  if (is.null(start_beta) || !is.numeric(start_beta)) {
    stop("start_beta must be a numeric vector")
  }
  if (any(is.na(y))) {
    stop("y must not contain missing values")
  }
  n <- length(y)

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
    if (is.null(names(start_delta))) {
      names(start_delta) <- paste0("delta", seq_along(start_delta))
    }
  } else {
    start_delta <- numeric(0)
  }

  start_params <- c(start_beta, start_delta)

  # Gauss-Hermite nodes/weights for N(0,1) integration (n = 20, hardcoded).
  z_nodes <- c(
    -7.6190485, -6.5105902, -5.5787388, -4.7345813, -3.9439674, -3.1890148,
    -2.4586636, -1.7452473, -1.0429453, -0.3469642, 0.3469642, 1.0429453,
    1.7452473, 2.4586636, 3.1890148, 3.9439674, 4.7345813, 5.5787388,
    6.5105902, 7.6190485
  )
  weights <- c(
    1.257801e-13, 2.482062e-10, 6.127490e-08, 4.402121e-06, 1.288263e-04,
    1.830103e-03, 1.399784e-02, 6.150637e-02, 1.617393e-01, 2.607931e-01,
    2.607931e-01, 1.617393e-01, 6.150637e-02, 1.399784e-02, 1.830103e-03,
    1.288263e-04, 4.402121e-06, 6.127490e-08, 2.482062e-10, 1.257801e-13
  )
  log_weights <- log(weights)

  X_w <- add_intercept(w)
  w_theta <- c(X_w %*% theta_hat)

  row_log_sum_exp <- function(mat) {
    row_max <- apply(mat, 1, max)
    row_max + log(rowSums(exp(mat - row_max)))
  }

  neg_loglik <- function(params) {
    beta <- params[seq_len(nb)]
    names(beta) <- beta_names
    if (length(params) > nb) {
      delta <- params[(nb + 1):length(params)]
      names(delta) <- names(start_delta)
    } else {
      delta <- numeric(0)
    }
    beta_g <- beta[["gf"]]
    beta_w <- beta[beta_names[-1]]

    base_linear <- beta_g * (posterior$a * gc + posterior$b * w_theta) +
      c(X_w %*% beta_w)

    # Accumulate log-likelihoods across quadrature nodes.
    loglik_mat <- matrix(NA_real_, nrow = n, ncol = length(z_nodes))
    for (j in seq_along(z_nodes)) {
      linpred <- base_linear + beta_g * posterior$c * z_nodes[j]
      loglik_vec <- loglik_fn(y = y, linpred = linpred, delta = delta)
      if (length(loglik_vec) != n) {
        stop("loglik_fn must return a vector of length length(y).")
      }
      if (any(!is.finite(loglik_vec))) {
        return(1e12)
      }
      loglik_mat[, j] <- loglik_vec + log_weights[j]
    }

    ll <- sum(row_log_sum_exp(loglik_mat))
    if (!is.finite(ll)) {
      return(1e12)
    }
    -ll
  }
  
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
  if (is.matrix(opt$hessian)) {
    vcov_try <- try(solve(opt$hessian), silent = TRUE)
    if (!inherits(vcov_try, "try-error")) {
      vcov_all <- vcov_try
      standard_errors <- sqrt(diag(vcov_all))[seq_len(nb)]
      names(standard_errors) <- beta_names
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
