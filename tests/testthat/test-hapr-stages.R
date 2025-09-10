test_that("calculate_beta is available and basic mapping works", {
  skip_if_not_installed("hapr")  # in case running in a different context

  expect_true(exists("calculate_beta", where = asNamespace("hapr"), inherits = FALSE))

  coeffs <- list(
    gamma = c("(Intercept)" = 0.1, gc = 0.2, w1 = 0.3),
    theta = c("(Intercept)" = 0.05, w1 = 0.25)
  )
  post <- list(a = 0.8, b = 0.2, c = 0.5)

  beta <- hapr:::calculate_beta("lm", coeffs, post)

  expect_true("gf" %in% names(beta))
  expect_false("gc" %in% names(beta))
  expect_equal(unname(beta["gf"]), 0.2 / 0.8, tolerance = 1e-12)
})

test_that("stage-1 -> stage-2 runs and is order-stable (lm)", {
  set.seed(123)

  n <- 300
  w <- data.frame(w1 = rnorm(n), w2 = rnorm(n))
  gc <- as.numeric(scale(rnorm(n)))
  y  <- 0.5 * gc + 0.3 * w$w1 - 0.2 * w$w2 + rnorm(n)

  fs <- hapr_first_stage(y = y, gc = gc, w = w, model_type = "lm")
  ss <- hapr_second_stage(fs, improvement_ratio = 1.2)

  expect_s3_class(fs, "hapr_first_stage_fit")
  expect_s3_class(ss, "hapr_fit")
  expect_true(is.numeric(ss$coefficients$beta))
  expect_true(is.numeric(ss$standard_errors))
  expect_true(is.matrix(ss$vcov_beta))
  expect_true(is.data.frame(ss$ci_beta))

  nb <- names(ss$coefficients$beta)
  ns <- names(ss$standard_errors)
  expect_setequal(nb, ns)

  z <- qnorm(0.975)
  est <- ss$ci_beta$Estimate
  se  <- ss$ci_beta$`Std.Error`
  low <- ss$ci_beta$Lower
  upp <- ss$ci_beta$Upper
  expect_equal(low, est - z*se, tolerance = 1e-7)
  expect_equal(upp, est + z*se, tolerance = 1e-7)

  # order stability
  w_perm <- w[, c("w2", "w1")]
  fs2 <- hapr_first_stage(y = y, gc = gc, w = w_perm, model_type = "lm")
  ss2 <- hapr_second_stage(fs2, improvement_ratio = 1.2)

  b1 <- ss$coefficients$beta[sort(names(ss$coefficients$beta))]
  b2 <- ss2$coefficients$beta[sort(names(ss2$coefficients$beta))]
  expect_equal(b2, b1, tolerance = 1e-7)

  se1 <- ss$standard_errors[sort(names(ss$standard_errors))]
  se2 <- ss2$standard_errors[sort(names(ss2$standard_errors))]
  expect_equal(se2, se1, tolerance = 1e-7)
})

test_that("bootstrap SDs are in the same ballpark as delta SEs (quick smoke test)", {
  set.seed(42)

  n <- 250
  w <- data.frame(w1 = rnorm(n), w2 = rnorm(n))
  gc <- as.numeric(scale(rnorm(n)))
  y  <- 0.6 * gc + 0.25 * w$w1 - 0.15 * w$w2 + rnorm(n)

  fs <- hapr_first_stage(y, gc, w, "lm")
  ss <- hapr_second_stage(fs, improvement_ratio = 1.2)

  B <- 100
  boot_betas <- replicate(B, {
    idx <- sample.int(n, replace = TRUE)
    fsb <- suppressWarnings(hapr_first_stage(y[idx], gc[idx], w[idx, ], "lm"))  # quiet numeric warning
    ssb <- hapr_second_stage(fsb, improvement_ratio = 1.2)
    ssb$coefficients$beta
  })

  emp_sd <- apply(boot_betas, 1, sd)
  del_se <- ss$standard_errors[names(emp_sd)]

  ratio <- emp_sd / del_se
  expect_true(all(is.finite(ratio)))
  expect_true(all(ratio > 0.5 & ratio < 2.0))
})
