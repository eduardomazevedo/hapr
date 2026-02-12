#' Heritability adjusted prediction
#'
#' @description
#' Fits the HARP model given the outcome y, PRS gc, and control variables w.
#' The improvement ratio specifies the expected increase in R-squared of the
#' future GWAS relative to the current GWAS.
#'
#' @param y Outcome variable
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#' @param model_type Type of model to fit ("lm" or "probit")
#' @param improvement_ratio The ratio of R-squared of the future fit to the current fit.
#' @param stage1_variance_method How to estimate stage-1 variance uncertainty for gc~w.
#'   One of "scaled_gc" (default) or "ols".
#'
#' @return A hapr_fit object.
#' @export
hapr <- function(y, gc, w, model_type, improvement_ratio, stage1_variance_method = c("scaled_gc", "ols")) {
  first_stage <- hapr_first_stage(y, gc, w, model_type, stage1_variance_method = stage1_variance_method)
  second_stage <- hapr_second_stage(first_stage, improvement_ratio)
  second_stage
}
