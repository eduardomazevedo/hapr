#' HAPR maximum likelihood estimation for survival models (analytic gradient)
#'
#' @description
#' Fits the HAPR parametric survival models by maximum likelihood using
#' Gauss-Hermite quadrature with analytic gradients. Supported models are
#' exponential and Weibull.
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
hapr_mle_survival_grad <- function(
    event_time,
    event_status,
    gc,
    w,
    improvement_ratio,
    model_type = "exponential",
    start_beta,
    start_delta = NULL,
    control = list()) {
  hapr_mle_survival(
    event_time = event_time,
    event_status = event_status,
    gc = gc,
    w = w,
    improvement_ratio = improvement_ratio,
    model_type = model_type,
    start_beta = start_beta,
    start_delta = start_delta,
    use_analytic_gradient = TRUE,
    control = control
  )
}
