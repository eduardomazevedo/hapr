#' Predict Outcomes from a hapr_fit Object
#'
#' @description
#' This function generates predictions using a fitted `hapr_fit` object.
#' It computes predictions based on specified covariates, supporting
#' linear (`lm`), `probit`, and `cox` regression.
#'
#' @param object An object of class `hapr_fit`, containing fitted regression models.
#' @param newdata A data frame with new observations for prediction.
#' @param covariates A character vector specifying which covariates to use for prediction.
#'        Defaults to `c('w', 'gc_w', 'gf_w')`: just covariates w, or combine them with gc or gf.
#' @param type A character string indicating the type of prediction.
#'        - For `lm`/`probit`, defaults to `"response"`.
#'        - For `cox`, can be `"lp"` (linear predictor) or `"response"`/`"risk"` (exp of linear predictor).
#' @param ... Additional arguments (currently ignored).
#' @return A data frame with the same rows as `newdata`, plus columns for predictions.
#'
#' @method predict hapr_fit
#' @export
predict.hapr_fit <- function(object, newdata, covariates = c("w", "gc_w", "gf_w"), type = "response", ...) {
  fit <- object
  # Validate model type
  if (!fit$model_type %in% c("lm", "probit", "cox")) {
    stop("Model type must be 'lm', 'probit', or 'cox'.")
  }

  # Prepare an output data frame
  results <- newdata

  # Loop over each requested covariate
  for (cov in covariates) {
    # Identify the corresponding model component by naming convention
    model_key <- paste0("y_on_", cov)

    # Skip if the regression for this covariate does not exist
    if (!model_key %in% names(fit$regressions)) {
      warning(sprintf("No regression found for covariate '%s'; skipping.", cov))
      next
    }

    # Extract regression object for the covariate
    reg_obj <- fit$regressions[[model_key]]
    
    # Extract coefficients and terms
    if (!"coefficients" %in% names(reg_obj)) {
      stop(sprintf("The regression for '%s' has no 'coefficients'.", cov))
    }
    coefs <- reg_obj$coefficients
    
    if (!"terms" %in% names(reg_obj)) {
      stop(sprintf("The regression for '%s' has no 'terms' object.", cov))
    }
    mt <- reg_obj$terms
    xlevels <- reg_obj$xlevels

    # Handle variable name mapping: if predicting with gf_w, terms object may reference gc
    # but coefficients and newdata use gf
    newdata_for_pred <- newdata
    if (cov == "gf_w" && "gf" %in% names(newdata) && !"gc" %in% names(newdata)) {
      # Temporarily rename gf to gc for model.frame (terms object expects gc)
      newdata_for_pred <- newdata
      newdata_for_pred$gc <- newdata_for_pred$gf
      newdata_for_pred$gf <- NULL
    }

    # Build model matrix using terms object (properly handles factor levels)
    # Delete response from terms for prediction
    mt_pred <- delete.response(mt)
    mf_new <- model.frame(mt_pred, data = newdata_for_pred, xlev = xlevels, na.action = na.pass)
    X <- model.matrix(mt_pred, data = mf_new)

    # Compute linear predictor - match coefficient names to column names
    # For gf_w, coefficients have "gf" but X has "gc", so map them
    X_colnames <- colnames(X)
    if (cov == "gf_w" && "gc" %in% X_colnames && "gf" %in% names(coefs)) {
      X_colnames[X_colnames == "gc"] <- "gf"
      colnames(X) <- X_colnames
    }
    
    common_vars <- intersect(names(coefs), colnames(X))
    X_sub <- X[, common_vars, drop = FALSE]
    coefs_sub <- coefs[common_vars]
    Xbeta <- as.vector(X_sub %*% coefs_sub)

    # Determine what to return based on model type and prediction type
    if (fit$model_type == "cox") {
      if (type == "lp") {
        results[[paste0("y_hat_", cov)]] <- Xbeta
      } else if (type %in% c("response", "risk")) {
        results[[paste0("y_hat_", cov)]] <- exp(Xbeta)
      } else {
        stop("For Cox models, `type` must be 'lp', 'response', or 'risk'.")
      }
    } else if (fit$model_type == "lm") {
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

  return(results)
}
