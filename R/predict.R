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

    # Extract the stripped model for the covariate
    model_obj <- fit$regressions[[model_key]]$stripped_model

    # Branch for Cox vs. LM/Probit
    if (fit$model_type == "cox") {
      # 1. Extract coefficients
      if (!"coefficients" %in% names(model_obj)) {
        stop(sprintf("The stripped Cox model for '%s' has no 'coefficients'.", cov))
      }
      coefs <- model_obj$coefficients

      # 2. Build model matrix (handles factors automatically)
      X <- model.matrix(~., data = newdata)

      # 3. Compute linear predictor from matching columns
      common_vars <- intersect(names(coefs), colnames(X))
      X_sub <- X[, common_vars, drop = FALSE]
      Xbeta <- as.vector(X_sub %*% coefs)

      # 4. Determine what to return
      if (type == "lp") {
        results[[paste0("y_hat_", cov)]] <- Xbeta
      } else if (type %in% c("response", "risk")) {
        results[[paste0("y_hat_", cov)]] <- exp(Xbeta)
      } else {
        stop("For Cox models, `type` must be 'lp', 'response', or 'risk'.")
      }
    } else {
      # Use standard `predict` for 'lm' or 'probit'
      results[[paste0("y_hat_", cov)]] <- predict(model_obj, newdata = newdata, type = type)
    }
  }

  return(results)
}
