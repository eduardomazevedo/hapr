#' Simulate data from a hapr_fit object
#'
#' @param object A hapr_fit object
#' @param w A data frame of covariates
#' @param gc Optional vector of genetic scores. If NULL, they will be simulated
#' @param repetitions Number of times to repeat each row in w
#'
#' @return A data frame with simulated values
#' @export
simulate.hapr_fit <- function(object, w, gc = NULL, repetitions = 1) {
  # Check if w is a data frame
  if (!is.data.frame(w)) {
    stop("The 'w' parameter must be a data frame.")
  }

  fit <- object
  n <- nrow(w) * repetitions
  w <- as.data.frame(w[rep(seq_len(nrow(w)), repetitions), ])

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

  # Simulate outcome based on model type
  y <- NULL
  if (fit$model_type == "lm") {
    beta <- fit$coefficients$beta
    data_for_simulation <- cbind(gf = gf, w)
    linear_predictor <- model.matrix(~., data = data_for_simulation) %*% beta
    
    # For linear models, add noise
    if (is.null(fit$additional_parameters$var_eta)) {
      # If var_eta not available, use a default noise level
      var_eta <- 1
    } else {
      var_eta <- fit$additional_parameters$var_eta
    }
    y <- linear_predictor + rnorm(n, 0, sqrt(var_eta))
    
  } else if (fit$model_type == "probit") {
    beta <- fit$coefficients$beta
    data_for_simulation <- cbind(gf = gf, w)
    linear_predictor <- model.matrix(~., data = data_for_simulation) %*% beta
    
    # For probit models, generate binary outcomes
    p <- pnorm(linear_predictor)
    y <- rbinom(n, 1, p)
    
  } else if (fit$model_type == "cox") {
    beta <- fit$coefficients$beta
    data_for_simulation <- cbind(gf = gf, w)
    linear_predictor <- model.matrix(~., data = data_for_simulation)[,-1] %*% beta[-1]
    
    # For cox models, simulate survival times - this is complex
    # and would require more detailed implementation based on baseline hazard
    # Here's a placeholder approach
    if (!is.null(fit$additional_parameters$baseline_hazard)) {
      # This is a simplified approach and might need to be expanded
      # based on specific requirements for Cox model simulation
      hazard_ratio <- exp(linear_predictor)
      # Simple exponential survival time simulation
      y <- rexp(n, hazard_ratio)
    } else {
      warning("Baseline hazard not available, Cox survival times not simulated")
      y <- rep(NA, n)
    }
  }

  # Return simulated data
  result <- data.frame(
    y = y,
    gf = gf, 
    gc = gc
  )
  
  # Add covariates
  result <- cbind(result, w)

  return(result)
}
