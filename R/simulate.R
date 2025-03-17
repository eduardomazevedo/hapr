# predict.R
simulate <- function(object, newdata, ...) {
  UseMethod("simulate")
}

#' Simulate Data from a hapr_fit Object
#'
#' This function simulates data based on a fitted `hapr_fit` model object.
#'
#' @param fit A `hapr_fit` object containing the model fit.
#' @param w A matrix of covariates.
#' @param gc A numeric vector of genetic components. If `NULL`, it will be generated internally.
#' @param repetitions An integer specifying the number of repetitions for the simulation. Default is 1.
#'
#' @return A data frame containing the simulated genetic factors (`gf`), genetic components (`gc`), and covariates (`w`).
#'
#' @examples
#' # Assuming `fit` is a hapr_fit object and `w` is a matrix of covariates
#' simulate.hapr_fit(fit, w, repetitions = 10)
#'
simulate.hapr_fit <- function(fit, w, gc = NULL, repetitions = 1) {
  n <- nrow(w) * repetitions
  w <- w[rep(seq_len(nrow(w)), repetitions), ]

  wtheta <- model.matrix(~., data = w) %*% fit$first_stage$gc_w_results$theta

  if (is.null(gc)) {
    gf <- rnorm(n, 0, sqrt(fit$second_stage$var_v)) + wtheta
    gc <- gf + rnorm(n, 0, sqrt(fit$second_stage$var_epsilon))
  } else {
    a <- fit$second_stage$posterior_parameters$a
    b <- fit$second_stage$posterior_parameters$b
    c <- fit$second_stage$posterior_parameters$c

    
    z <- rnorm(n, 0, 1)
    gf <- a * gc + b * wtheta + c * z
  }

  # Bind gf and gc as new columns to the data frame w
  result <- cbind(gf = gf, gc = gc, w)

  return(result)
}