#' HAPR probit model first stage fit
#' @param y Outcome variable
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#'
#' @return A hapr_probit_first_stage_fit object containing the results of the first stage.
#' @details
#' Fits the HARP model given the outcome y, PRS gc, and control variables w. This returns
#' a first stage fit, which does not need to assume an improvement ratio. Run
#' hapr_probit_second_stage(first_stage_fit, improvement_ratio) to specify an improvement ratio and get the full model.
#' @export
hapr_probit_first_stage <- function(y, gc, w) {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = "probit")
  y <- preprocessed$y
  gc <- preprocessed$gc
  w <- preprocessed$w
  rm(preprocessed)

  # Get regression results of gc on w
  gc_w_results <- gc_regression(gc, w)

  # Regress y on gc and w
  y_gc_w_results <- feasible_regression_probit(y, gc, w)

  # Compute max_r2_gf
  stats_list <- list(
    max_r2_gf = gc_w_results$max_improvement_ratio * y_gc_w_results$r2_gc
  )

  # Return
  result <- list(
    gc_w_results = gc_w_results,
    y_gc_w_results = y_gc_w_results,
    first_stage_stats = stats_list
  )
  class(result) <- "hapr_probit_first_stage_fit"
  result
}
