simulate <- function(object, w, gc = NULL, repetitions = 1) {
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
