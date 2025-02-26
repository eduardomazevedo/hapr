#' Heritability adjusted prediction for the Cox proportional hazards model
#'
#' Fits the HARP model given the outcome y, PRS gc, and control variables w.
#'
#' @param y Outcome variable: a Surv object compatible with survival::coxph()
#' @param gc Polygenic risk score (has to be normalized)
#' @param w Control variables
#' @param improvement_ratio The ratio of R-squared of the future fit to the current fit.
#' @param r2_current The R-squared of the current fit
#' @param r2_future The R-squared of the future fit
#'
#' @return A hapr_cox_fit object.
#' @export
hapr_cox <- function(y, gc, w, improvement_ratio = NULL, r2_current = NULL, r2_future = NULL) {
  first_stage <- hapr_cox_first_stage(y, gc, w)
  second_stage <- hapr_cox_second_stage(first_stage, improvement_ratio, r2_current, r2_future)
  second_stage
}
