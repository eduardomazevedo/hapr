#' Simulate Data from a hapr_fit Object
#'
#' @param object A `hapr_fit` object containing the model fit.
#' @param w A matrix or data frame of covariates.
#' @param gc A numeric vector of genetic components. If `NULL`, generated internally.
#' @param repetitions Number of repetitions for the simulation. Default is 1.
#' @param seed Random seed for reproducibility. Default is NULL.
#'
#' @return A data frame containing simulated genetic factors (`gf`), genetic components (`gc`), and covariates (`w`).
#'
#' @examples
#' simulate(fit, w, repetitions = 10, seed = 123)
#' @export
simulate <- function(object, ...) {
  UseMethod("simulate")
}

#' @export
simulate.hapr_fit <- function(object, w, gc = NULL, repetitions = 1, seed = NULL, ...) {
  if (!inherits(object, "hapr_fit")) stop("object must be of class 'hapr_fit'")
  if (!is.matrix(w) && !is.data.frame(w)) stop("w must be a matrix or data frame")

  if (!is.null(seed)) set.seed(seed)

  fit <- object
  n <- nrow(w) * repetitions
  w <- as.data.frame(w)[rep(seq_len(nrow(w)), repetitions), ]

  wtheta <- model.matrix(~., data = w) %*% fit$coefficients$theta

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
