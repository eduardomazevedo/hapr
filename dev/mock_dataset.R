#' Simulate data for linear-model HAPR examples
#'
#' @param n Number of observations
#' @param var_v Variance of latent v component
#' @param var_epsilon Variance of PRS measurement error
#' @param beta_g True coefficient on gf
#' @param beta_w Coefficients for intercept and w covariates
#' @param theta Coefficients for gc ~ w regression (including intercept)
#' @param var_y Variance of the outcome error term
#'
#' @return List with w, gf, gc, y
mock_dataset_lm <- function(
    n,
    var_v,
    var_epsilon,
    beta_g,
    beta_w,
    theta,
    var_y
) {
    p <- length(theta) - 1
    w = matrix(rnorm(n * p), n, p)
    e = rnorm(n, 0, sqrt(var_epsilon))
    v = rnorm(n, 0, sqrt(var_v))

    # Normalize w so that var(theta %*% w) = 1 - var_v - var_epsilon
    # Since w is standard Gaussian, var(theta'w) = theta'*theta = sum(theta^2)
    w = w * sqrt((1 - var_v - var_epsilon) / sum(theta[-1]^2))
    
    # Create w with intercept for internal calculations
    w_with_int <- cbind(1, w)
    
    gf = v + w_with_int %*% theta
    gc = gf + e
    y = beta_g * gf + w_with_int %*% beta_w + rnorm(n, 0, var_y)

    # Return w without intercept (as expected by hapr_first_stage)
    return(list(w = w, gf = gf, gc = gc, y = y))
}


#' Simulate data for probit-model HAPR examples
#'
#' @inheritParams mock_dataset_lm
#'
#' @return List with w, gf, gc, y (logical)
mock_dataset_probit <- function(
    n,
    var_v,
    var_epsilon,
    beta_g,
    beta_w,
    theta
) {
    results <- mock_dataset_lm(n, var_v, var_epsilon, beta_g, beta_w, theta, 1)
    y = results$y > 0
    return(list(w = results$w, gf = results$gf, gc = results$gc, y = y))
}

#' Simulate data for exponential survival HAPR examples
#'
#' @inheritParams mock_dataset_lm
#' @param censor_rate Rate parameter for exponential censoring (0 for none)
#'
#' @return List with w, gf, gc, event_time, event_status
mock_dataset_survival_exponential <- function(
    n,
    var_v,
    var_epsilon,
    beta_g,
    beta_w,
    theta,
    censor_rate = 0.0
) {
    base <- mock_dataset_lm(n, var_v, var_epsilon, beta_g, beta_w, theta, 0)
    w_with_int <- cbind(1, base$w)
    linpred = beta_g * base$gf + w_with_int %*% beta_w
    rate = exp(linpred)
    event_time = rexp(n, rate = rate)

    if (censor_rate > 0) {
        censor_time = rexp(n, rate = censor_rate)
        event_status = as.numeric(event_time <= censor_time)
        event_time = pmin(event_time, censor_time)
    } else {
        event_status = rep(1, n)
    }

    return(list(
        w = base$w,
        gf = base$gf,
        gc = base$gc,
        event_time = event_time,
        event_status = event_status
    ))
}
