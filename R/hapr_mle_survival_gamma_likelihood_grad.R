make_hapr_mle_likelihood_survival_gamma_grad <- function(
    event_time,
    event_status,
    gc,
    X_w,
    posterior,
    model_type,
    use_openmp = TRUE) {
  event_idx <- which(event_status == 1)
  censor_idx <- which(event_status == 0)

  event_time_event <- event_time[event_idx]
  gc_event <- gc[event_idx]
  X_w_event <- X_w[event_idx, , drop = FALSE]

  censor_time <- event_time[censor_idx]
  gc_censor <- gc[censor_idx]
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

  post_c_over_a <- posterior$c / posterior$a

  nll_grad <- function(params) {
    hapr_mle_survival_gamma_nll_split_grad_cpp(
      params = params,
      event_time = event_time_event,
      gc_event = gc_event,
      X_w_event = X_w_event,
      censor_time = censor_time,
      gc_censor = gc_censor,
      X_w_censor = X_w_censor,
      post_c_over_a = post_c_over_a,
      model_type = model_id,
      use_openmp = use_openmp
    )
  }

  list(
    fn = function(params) nll_grad(params)$value,
    gr = function(params) nll_grad(params)$gradient
  )
}
