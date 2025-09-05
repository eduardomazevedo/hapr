hapr_sim_table_boot <- function(n_sims = 200,
                                n      = 2000,
                                improvement_ratio = 1.5,
                                B_boot = 200,
                                seed   = 12345) {
  stopifnot(is.function(hapr_first_stage), is.function(hapr_second_stage), is.function(abc))
  set.seed(seed)

  # choose true beta (future coefficients)
  w_names <- paste0("w", 1:5)
  beta_true <- c("(Intercept)" = 0,
                 "gf" = 0.42,
                 setNames(c(0.10, -0.06, 0.04, 0.03, -0.02), w_names))
  terms_w <- setdiff(names(beta_true), c("(Intercept)","gf"))

  ## improvement ratio -> (a,b,c)
  var_epsilon <- 1 - 1/improvement_ratio
  var_v_plus_var_epsilon <- 0.50
  var_v <- var_v_plus_var_epsilon - var_epsilon
  if (var_v <= 0) stop("Choose a larger var_v_plus_var_epsilon or a smaller improvement_ratio.")
  post <- abc(var_epsilon = var_epsilon, var_v = var_v)
  a <- post$a; b <- post$b

  ## theta base (scaled to keep Var(gc) ~ 1)
  theta_base <- setNames(seq(0.15, by = -0.03, length.out = length(terms_w)), terms_w)

  scale_theta <- function(theta0, target_var) {
    s2 <- sum(theta0^2); if (s2 == 0) return(theta0)
    theta0 * sqrt(target_var / s2)
  }

  estimates_list <- vector("list", n_sims)
  se_boot_list   <- vector("list", n_sims)  # bootstrap SEs
  covers_list    <- vector("list", n_sims)
  terms_order    <- NULL
  true_beta_vec  <- NULL

  for (s in seq_len(n_sims)) {
    ## --- simulate one dataset from the *first-stage* model y ~ gc + w ---
    w <- as.data.frame(matrix(rnorm(n * length(terms_w)), n, length(terms_w)))
    names(w) <- terms_w

    theta_scaled <- scale_theta(theta_base, target_var = 1 - var_v_plus_var_epsilon)
    u  <- rnorm(n, sd = sqrt(var_v_plus_var_epsilon))
    gc <- as.numeric(as.matrix(w) %*% theta_scaled + u)

    gamma_gc <- a * beta_true["gf"]
    gamma_w  <- beta_true[terms_w] + as.numeric(beta_true["gf"]) * b * theta_scaled[terms_w]
    y <- 0 + gamma_gc * gc + as.numeric(as.matrix(w) %*% gamma_w) + rnorm(n, sd = 1)

    # fit once
    fs <- hapr_first_stage(y = y, gc = gc, w = w, model_type = "lm")
    ss <- hapr_second_stage(fs, improvement_ratio = improvement_ratio)

    beta_hat <- ss$coefficients$beta

    ## fix term order on first iteration
    if (is.null(terms_order)) {
      terms_order <- names(beta_hat)
      true_beta_vec <- beta_true
      missing <- setdiff(terms_order, names(true_beta_vec))
      if (length(missing)) {
        true_beta_vec <- c(true_beta_vec, stats::setNames(rep(0, length(missing)), missing))
      }
      true_beta_vec <- true_beta_vec[terms_order]
    }

    # bootstrap inside this dataset to get SE and CI for coverage
    boot_betas <- replicate(B_boot, {
      idx <- sample.int(n, replace = TRUE)
      fsb <- hapr_first_stage(y[idx], gc[idx], w[idx, , drop = FALSE], "lm")
      ssb <- hapr_second_stage(fsb, improvement_ratio = improvement_ratio)
      ssb$coefficients$beta[terms_order]
    })

    se_boot <- apply(boot_betas, 1, sd)
    z <- qnorm(0.975)
    L <- beta_hat[terms_order] - z * se_boot
    U <- beta_hat[terms_order] + z * se_boot
    cov_vec <- as.numeric(L <= true_beta_vec & true_beta_vec <= U)

    estimates_list[[s]] <- beta_hat[terms_order]
    se_boot_list[[s]]   <- se_boot[terms_order]
    covers_list[[s]]    <- cov_vec
  }

  estimates_mat <- do.call(rbind, estimates_list)
  seboot_mat    <- do.call(rbind, se_boot_list)
  covers_mat    <- do.call(rbind, covers_list)

  data.frame(
    parameter           = terms_order,
    true_beta           = as.numeric(true_beta_vec),
    mean_beta_hat       = colMeans(estimates_mat),
    stdev_beta_hat      = apply(estimates_mat, 2, sd),
    mean_standard_error = colMeans(seboot_mat),
    coverage_percentage = 100 * colMeans(covers_mat),
    row.names = NULL,
    check.names = FALSE
  )
}

# print table
library(hapr); devtools::load_all()
tbl <- hapr_sim_table_boot(n_sims = 100, n = 3000, improvement_ratio = 1.5, B_boot = 200)
print(tbl, row.names = FALSE)

