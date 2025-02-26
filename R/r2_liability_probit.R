#' Calculate liability-scale R-squared for probit models
#' 
#' Computes the liability-scale R-squared for a probit regression model. This measure 
#' represents the proportion of variance explained on the liability (latent) scale.
#'
#' @noRd
#' @param model A fitted model object of class "glm" with probit link function
#'
#' @return A numeric value between 0 and 1 representing the liability-scale R-squared.
#'   Returns 0 if the model only contains an intercept.
#'
#' @details 
#' The liability-scale R-squared is calculated as:
#' var(X\beta) / (var(X\beta) + 1)
#' where X\beta is the linear predictor excluding the intercept.
#'
#' @examples
#' \dontrun{
#' # Fit a probit model
#' model <- glm(y ~ x, family = binomial(link = "probit"))
#' r2 <- r2_liability_probit(model)
#' }
r2_liability_probit <- function(model) {
  # Check if the model is a glm object
  if (!inherits(model, "glm")) {
    stop("Error: The input model is not a glm object.")
  }
  
  # Check if the model uses a probit link function
  if (model$family$link != "probit") {
    stop("Error: The model does not use a probit link function.")
  }
  
  # Extract coefficients and design matrix (excluding intercept)
  X <- model.matrix(model)[, -1, drop = FALSE]  # Remove intercept
  beta_hat <- coef(model)[-1]                   # Remove intercept
  
  # Check if there are predictors (if not, return R^2 = 0)
  if (length(beta_hat) == 0) {
    warning("The model only has an intercept. R^2_liability is 0.")
    return(0)
  }
  
  # Compute variance of the linear predictor
  var_linear_predictor <- var(X %*% beta_hat)
  
  # Compute liability-scale R^2
  r2_liability <- var_linear_predictor / (var_linear_predictor + 1)
  
  return(r2_liability)
}
