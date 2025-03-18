#' Predict Outcomes from a hapr_fit Object
#'
#' This function generates predictions using a fitted `hapr_fit` object.
#' It computes predictions based on specified covariates.
#'
#' @param fit An object of class `hapr_fit`, containing fitted regression models.
#' @param newdata A data frame with new observations for prediction.
#' @param covariates A character vector specifying which covariates to use for prediction.
#'        Defaults to `c('w', 'gc_w', 'gf_w')`.
#' @param type A character string indicating the type of prediction. Defaults to "response".
#'
#' @return A data frame with the same structure as `newdata`, with additional columns for predicted values.
#'
#' @export
predict.hapr_fit <- function(fit, newdata, covariates = c('w', 'gc_w', 'gf_w'), type = "response") {
    # Check that fit$model_type is either "lm" or "probit"
    if (!fit$model_type %in% c("lm", "probit")) {
        stop("Model type must be either 'lm' or 'probit'.")
    }
    results <- newdata
    if ('w' %in% covariates) {
        results$y_hat_w <- predict(fit$regressions$y_on_w$stripped_model, newdata = newdata, type = type)
    }

    if ('gc_w' %in% covariates) {
        results$y_hat_gc_w <- predict(fit$regressions$y_on_gc_w$stripped_model, newdata = newdata, type = type)
    }

    if ('gf_w' %in% covariates) {
        results$y_hat_gf_w <- predict(fit$regressions$y_on_gf_w$stripped_model, newdata = newdata, type = type)
    }

    return(results)
}