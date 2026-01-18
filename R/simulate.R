#' Simulate data from a hapr_fit model
#'
#' @description
#' This function simulates data from a fitted HAPR model. It generates genetic scores
#' and outcomes based on the model parameters estimated in the hapr_fit object.
#'
#' @param object A hapr_fit object from hapr() or hapr_second_stage()
#' @param w A data frame of control variables/covariates
#' @param gc Optional vector of polygenic risk scores. If NULL, they will be simulated
#'        based on the parameters in the model
#' @param repetitions Number of times to repeat each row in w
#'
#' @details
#' The function uses the parameters from the HAPR model to simulate genetic data and outcomes.
#' When gc is not provided, both the true genetic effect (gf) and measured genetic score (gc)
#' are simulated according to the heritability parameters in the model. If gc is provided,
#' the true genetic effect (gf) is simulated.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item gf: Simulated true genetic value
#'   \item gc: Simulated or provided polygenic risk score
#'   \item Additional columns from the w data frame
#' }
#'
#' @examples
#' \dontrun{
#' # Fit a HAPR model
#' fit <- hapr(y = outcome, gc = genetic_score, w = covariates, 
#'             model_type = "lm", improvement_ratio = 2)
#'             
#' # Simulate 100 observations using the same covariates
#' sim_data <- hapr_simulate(fit, w = covariates, repetitions = 10)
#' 
#' # Simulate with provided genetic scores
#' sim_data2 <- hapr_simulate(fit, w = covariates, gc = genetic_score)
#' }
#'
#' @seealso \code{\link{hapr}}, \code{\link{hapr_first_stage}}, \code{\link{hapr_second_stage}}
#' @export
hapr_simulate <- function(object, w, gc = NULL, repetitions = 1) {
  # Check if w is a data frame
  if (!is.data.frame(w)) {
    stop("The 'w' parameter must be a data frame.")
  }
  # Check that object is a hapr_fit
  if (!inherits(object, "hapr_fit")) {
    stop("The 'object' parameter must be a hapr_fit object.")
  }
  # If gc is provided and repetitions is greater than 1, throw an error
  if (!is.null(gc) && repetitions > 1) {
    stop("If gc is provided, repetitions must be 1.")
  }

  fit <- object
  n <- nrow(w) * repetitions
  w <- as.data.frame(w[rep(seq_len(nrow(w)), repetitions), ])

  wtheta <- model.matrix(~., data = w) %*% fit$parameters$theta

  if (is.null(gc)) {
    gf <- rnorm(n, 0, sqrt(fit$stats$var_v)) + wtheta
    gc <- gf + rnorm(n, 0, sqrt(fit$stats$var_epsilon))
  } else {
    a <- fit$stats$posterior$a
    b <- fit$stats$posterior$b
    c <- fit$stats$posterior$c
    z <- rnorm(n, 0, 1)
    gf <- a * gc + b * wtheta + c * z
  }

  result <- cbind(gf = gf, gc = gc, w)

  return(result)
}
