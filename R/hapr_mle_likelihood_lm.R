#' Build linear-model MLE objective for HAPR
#'
#' @param y Outcome vector
#' @param gc Normalized PRS vector
#' @param w_theta Precomputed w %*% theta vector
#' @param X_w Design matrix with intercept and w columns
#' @param posterior List with a, b, c parameters
#'
#' @return Function(params) returning negative log-likelihood
make_hapr_mle_likelihood_lm <- function(y, gc, w_theta, X_w, posterior) {
  avg_linpred <- posterior$a * gc + posterior$b * w_theta
  function(params) {
    hapr_mle_lm_nll_cpp(
      params = params,
      y = y,
      avg_linpred = avg_linpred,
      X_w = X_w,
      post_c = posterior$c
    )
  }
}
