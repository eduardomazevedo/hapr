test_that("abc calculates correct values", {
  # Test with known values
  result <- abc(var_epsilon = 4, var_v = 9)

  expect_equal(result$a, (1 / 4) / (1 / 4 + 1 / 9)) # a = precision_epsilon / (precision_epsilon + precision_v)
  expect_equal(result$b, (1 / 9) / (1 / 4 + 1 / 9)) # b = precision_v / (precision_epsilon + precision_v)
  expect_equal(result$c, 1 / sqrt(1 / 4 + 1 / 9)) # c = 1 / sqrt(precision_sum)

  # Edge case: var_epsilon and var_v are equal
  result_equal <- abc(var_epsilon = 5, var_v = 5)
  expect_equal(result_equal$a, 0.5)
  expect_equal(result_equal$b, 0.5)
  expect_equal(result_equal$c, 1 / sqrt(2 / 5))
})
