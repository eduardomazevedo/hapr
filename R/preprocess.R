#' Preprocess inputs for HAPR models
#'
#' @param y Outcome variable. For linear models, must be numeric. For probit models,
#'   must be numeric, logical, or factor with exactly 2 levels. For cox models,
#'   must be a Surv object compatible with survival::coxph().
#' @param gc Polygenic risk score, will be normalized
#' @param w Data frame of control variables
#' @param model_type Type of model to fit. Must be one of: "lm" (linear model),
#'   "probit" (probit model), or "cox" (cox proportional hazards model)
#'
#' @return A list containing the preprocessed variables:
#'   \item{y}{Processed outcome variable}
#'   \item{gc}{Normalized polygenic risk score}
#'   \item{w}{Processed control variables}
#'
#' @details
#' The function:
#' 1. Validates inputs based on model type
#' 2. Normalizes the polygenic risk score
#' 3. Handles missing values by removing incomplete cases
#' 4. Ensures consistent dimensions across inputs
#'
#' @noRd
preprocess <- function(y, gc, w, model_type) {
  # Check that model_type is valid
  if (!model_type %in% c("lm", "probit", "cox")) {
    stop("model_type must be one of: 'lm', 'probit', 'cox'")
  }
  # Process y
  if (model_type == "lm") {
    y <- as.numeric(y)
  } else if (model_type == "probit") {
    if (!is.numeric(y) && !is.logical(y) && !is.factor(y)) {
      stop("For probit models, y must be numeric, logical, or factor")
    }
    if (is.numeric(y) || is.logical(y)) {
      y <- as.factor(y)
    }
    if (nlevels(y) != 2) {
      stop("For probit models, y must have exactly 2 levels")
    }
  } else if (model_type == "cox") {
    if (!survival::is.Surv(y)) {
      stop("For cox models, y must be a Surv object. Use survival::Surv() to create one.")
    }

    cox_test <- try(survival::coxph(y ~ 1), silent = TRUE)
    if (inherits(cox_test, "try-error")) {
      stop("For cox models, y must be a Surv object compatible with survival::coxph().")
    }
  }

  # Ensure no N/A values are being supplied to the HAPR function.
  if (any(is.na(y))) {
    stop("There are N/A values in the outcome variable you are sending to HAPR.")
  }
  if (any(is.na(gc))) {
    stop("There are N/A values in the polygenic scores you are sending to HAPR.")
  }
  if (any(is.na(w))) {
    stop("There are N/A values in the covariates you are sending to HAPR.")
  }
  

  # Check that y, gc, and w have the same number of observations
  if (length(y) != length(gc) || length(y) != nrow(w)) {
    stop("y, gc, and w must have the same number of observations")
  }

  # Scale gc and make sure it is numeric
  gc <- gc |>
    as.numeric() |>
    scale() |>
    as.numeric()

  # Ensure w is a data frame
  if (!is.data.frame(w)) {
    stop("w must be a data frame.")
  }

  # Make sure w has no variables called gc or gf
  if (any(names(w) %in% c("gc", "gf"))) {
    stop("w must not have variables called gc or gf")
  }


  list(y = y, gc = gc, w = w)
}
