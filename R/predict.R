#' Predict Outcomes from a hapr_fit Object
#'
#' @description
#' This function generates predictions using a fitted `hapr_fit` object.
#' It computes predictions based on specified covariates, supporting
#' linear (`lm`) and `probit` regression. Predictions require numeric
#' design matrices; no formula-based preprocessing is performed.
#'
#' @param object An object of class `hapr_fit`, containing fitted regression models.
#' @param newdata Either a numeric matrix of covariates `w`, or a list with
#'   elements `w` (numeric matrix), and optionally `gc` and `gf` vectors. The
#'   columns of `w` are renamed to `w1`, `w2`, ... to match `hapr()` behavior.
#' @param covariates A character vector specifying which covariates to use for prediction.
#'        Defaults to `c('w', 'gc_w', 'gf_w')`: just covariates w, or combine them with gc or gf.
#' @param type A character string indicating the type of prediction.
#'        - For `lm`, defaults to `"response"`.
#'        - For `probit`, can be `"response"` or `"link"` (linear predictor).
#' @param ... Additional arguments (currently ignored).
#' @return A data frame with the same rows as `newdata`, plus columns for predictions.
#'
#' @method predict hapr_fit
#' @export
predict.hapr_fit <- function(object, newdata, covariates = c("w", "gc_w", "gf_w"), type = "response", ...) {
  fit <- object
  # Validate model type
  if (!fit$model_type %in% c("lm", "probit")) {
    stop("Model type must be 'lm' or 'probit'.")
  }

  if (is.matrix(newdata)) {
    w <- newdata
    gc <- NULL
    gf <- NULL
    results <- data.frame(w)
  } else if (is.list(newdata)) {
    w <- newdata$w
    gc <- newdata$gc
    gf <- newdata$gf
    results <- as.data.frame(newdata)
  } else {
    stop("newdata must be a numeric matrix (w) or a list with w/gc/gf.")
  }

  if (!is.matrix(w) || !is.numeric(w)) {
    stop("newdata$w must be a numeric matrix.")
  }
  colnames(w) <- paste0("w", seq_len(ncol(w)))
  X_w <- cbind(1, w)
  colnames(X_w)[1] <- "(Intercept)"

  # Loop over each requested covariate
  for (cov in covariates) {
    if (cov == "w") {
      if (is.null(fit$regressions$y_on_w$coefficients)) {
        warning("No regression found for covariate 'w'; skipping.")
        next
      }
      coefs <- fit$regressions$y_on_w$coefficients
      if (!all(names(coefs) %in% colnames(X_w))) {
        stop("Covariate names for 'w' do not match coefficient names.")
      }
      Xbeta <- as.vector(X_w[, names(coefs), drop = FALSE] %*% coefs)
    } else if (cov == "gc_w") {
      if (is.null(gc)) {
        stop("newdata must include gc for covariates = 'gc_w'.")
      }
      if (is.null(fit$regressions$y_on_gc_w$coefficients)) {
        warning("No regression found for covariate 'gc_w'; skipping.")
        next
      }
      coefs <- fit$regressions$y_on_gc_w$coefficients
      X_gc_w <- cbind(gc = gc, X_w)
      if (!all(names(coefs) %in% colnames(X_gc_w))) {
        stop("Covariate names for 'gc_w' do not match coefficient names.")
      }
      Xbeta <- as.vector(X_gc_w[, names(coefs), drop = FALSE] %*% coefs)
    } else if (cov == "gf_w") {
      if (is.null(gf)) {
        stop("newdata must include gf for covariates = 'gf_w'.")
      }
      coefs <- fit$parameters$beta
      X_gf_w <- cbind(gf = gf, X_w)
      if (!all(names(coefs) %in% colnames(X_gf_w))) {
        stop("Covariate names for 'gf_w' do not match coefficient names.")
      }
      Xbeta <- as.vector(X_gf_w[, names(coefs), drop = FALSE] %*% coefs)
    } else {
      warning(sprintf("Unknown covariate '%s'; skipping.", cov))
      next
    }

    # Determine what to return based on model type and prediction type
    if (fit$model_type == "lm") {
      if (type != "response") {
        stop("For lm models, `type` must be 'response'.")
      }
      results[[paste0("y_hat_", cov)]] <- Xbeta
    } else if (fit$model_type == "probit") {
      if (!type %in% c("response", "link")) {
        stop("For probit models, `type` must be 'response' or 'link'.")
      }
      if (type == "link") {
        results[[paste0("y_hat_", cov)]] <- Xbeta
      } else if (type == "response") {
        results[[paste0("y_hat_", cov)]] <- pnorm(Xbeta)
      }
    }
  }

  results
}
