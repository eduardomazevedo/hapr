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
    w <- cbind(1, w)

    gf = v + w %*% theta
    gc = gf + e
    y = beta_g * gf + w %*% beta_w + rnorm(n, 0, var_y)

    return(list(w = w, gf = gf, gc = gc, y = y))
}


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