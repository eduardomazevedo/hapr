#' HAPR Cox proportional hazards model first stage fit
#' @param y Outcome variable: a Surv object compatible with survival::coxph()
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#'
#' @return A hapr_cox_first_stage_fit object containing the results of the first stage.
#' @details
#' Fits the HARP model given the outcome y, PRS gc, and control variables w. This returns
#' a first stage fit, which does not need to assume an improvement ratio. Run
#' hapr_cox_second_stage(first_stage_fit, improvement_ratio) to specify an improvement ratio and get the full model.
#' @export
hapr_cox_first_stage <- function(y, gc, w) {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = "cox")
  y <- preprocessed$y
  gc <- preprocessed$gc
  w <- preprocessed$w
  rm(preprocessed)

  # Get regression results of gc on w
  gc_w_results <- gc_regression(gc, w)

  # Regress y on gc and w
  y_gc_w_results <- feasible_regression_cox(y, gc, w)

  # Return
  result <- c(gc_w_results, y_gc_w_results)
  class(result) <- "hapr_cox_first_stage_fit"
  result
}
