#' HAPR linear model first stage fit
#' @param y Outcome variable
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#'
#' @return A hapr_lm_first_stage_fit object containing the results of the first stage.
#' @details
#' Fits the HARP model given the outcome y, PRS gc, and control variables w. This returns
#' a first stage fit, which does not need to assume an improvement ratio. Run
#' extrapolate(first_stage_fit) to specify an improvement ratio and get the full model.
#' @export
hapr_lm_first_stage <- function(y, gc, w) {
  # Preprocess inputs
  preprocessed <- preprocess(y, gc, w, model_type = "lm")
  y <- preprocessed$y
  gc <- preprocessed$gc
  w <- preprocessed$w
  rm(preprocessed)

  # Get regression results of gc on w
  gc_w_results <- regress_gc_on_w(gc, w)

  # Regress y on gc and w
  y_gc_w_results <- regress_lm_y_on_gc_and_w(y, gc, w)

  # Compute max_r2_gf
  max_r2_gf_list <- list(
    max_r2_gf = gc_w_results$max_improvement_ratio * y_gc_w_results$r2_gc
  )

  # Return
  result <- c(gc_w_results, y_gc_w_results, max_r2_gf_list)
  class(result) <- "hapr_lm_first_stage_fit"
  return(result)
}