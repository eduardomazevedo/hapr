#' Calculate a, b, and c
#' Calculates the constants of the posterior distribution of gf.
#' @noRd
#' @param var_epsilon Variance of epsilon
#' @param var_v Variance of V
#' @return A list containing a, b, and c
abc <- function(var_epsilon, var_v) {
  # Calculate precision
  precision_epsilon <- 1 / var_epsilon
  precision_v <- 1 / var_v

  # Calculate denominator (1/c^2)
  precision_sum <- precision_epsilon + precision_v

  # Calculate a and b
  a <- precision_epsilon / precision_sum
  b <- precision_v / precision_sum

  # Calculate c
  c <- 1 / sqrt(precision_sum)

  # Return results
  list(
    a = a,
    b = b,
    c = c
  )
}
