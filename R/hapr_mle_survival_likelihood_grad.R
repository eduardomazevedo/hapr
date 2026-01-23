make_hapr_mle_likelihood_survival_grad <- function(
    event_time,
    event_status,
    gc,
    w_theta,
    X_w,
    posterior,
    model_type,
    use_openmp = TRUE,
    use_optimized_gradient = FALSE) {
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
  if (!is.logical(use_openmp) || length(use_openmp) != 1) {
    stop("use_openmp must be a single logical value.")
  }
  if (!is.logical(use_optimized_gradient) || length(use_optimized_gradient) != 1) {
    stop("use_optimized_gradient must be a single logical value.")
  }

  nll_grad <- function(params) {
    if (use_optimized_gradient) {
      hapr_mle_survival_nll_split_grad_cpp_optimized(
        params = params,
        event_time = event_time_event,
        avg_linpred_event = avg_event,
        X_w_event = X_w_event,
        censor_time = censor_time,
        avg_linpred_censor = avg_censor,
        X_w_censor = X_w_censor,
        post_c = posterior$c,
        model_type = model_id,
        use_openmp = use_openmp
      )
    } else {
      hapr_mle_survival_nll_split_grad_cpp(
        params = params,
        event_time = event_time_event,
        avg_linpred_event = avg_event,
        X_w_event = X_w_event,
        censor_time = censor_time,
        avg_linpred_censor = avg_censor,
        X_w_censor = X_w_censor,
        post_c = posterior$c,
        model_type = model_id,
        use_openmp = use_openmp
      )
    }
  }

  list(
    fn = function(params) nll_grad(params)$value,
    gr = function(params) nll_grad(params)$gradient
  )
}
