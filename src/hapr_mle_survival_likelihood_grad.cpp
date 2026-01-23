#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
constexpr int kNumNodes = 20;
static const double kNodes[kNumNodes] = {
  -7.6190485, -6.5105902, -5.5787388, -4.7345813, -3.9439674, -3.1890148,
  -2.4586636, -1.7452473, -1.0429453, -0.3469642, 0.3469642, 1.0429453,
  1.7452473, 2.4586636, 3.1890148, 3.9439674, 4.7345813, 5.5787388,
  6.5105902, 7.6190485
};
static const double kLogWeights[kNumNodes] = {
  -29.704241251, -22.116761264, -16.607895540, -12.333424088, -8.957045573,
  -6.303383030, -4.268852247, -2.788614532, -1.821769499, -1.344027906,
  -1.344027906, -1.821769499, -2.788614532, -4.268852247, -6.303383030,
  -8.957045573, -12.333424088, -16.607895540, -22.116761264, -29.704241251
};
}  // namespace

// [[Rcpp::export]]
Rcpp::List hapr_mle_survival_nll_split_grad_cpp(
    const arma::vec& params,
    const arma::vec& event_time,
    const arma::vec& avg_linpred_event,
    const arma::mat& X_w_event,
    const arma::vec& censor_time,
    const arma::vec& avg_linpred_censor,
    const arma::mat& X_w_censor,
    double post_c,
    int model_type,
    bool use_openmp) {
  auto log_exp_density = [](double time, double linpred) {
    const double rate = std::exp(linpred);
    return linpred - rate * time;
  };
  auto log_exp_tail_probability = [](double time, double linpred) {
    const double rate = std::exp(linpred);
    return -rate * time;
  };
  auto dlog_exp_density_dlinpred = [](double time, double linpred) {
    const double rate = std::exp(linpred);
    return 1.0 - rate * time;
  };
  auto dlog_exp_tail_dlinpred = [](double time, double linpred) {
    const double rate = std::exp(linpred);
    return -rate * time;
  };
  auto log_weibull_density = [](double time, double linpred, double log_k, double k) {
    const double log_time = std::log(time);
    const double log_scaled_time = log_time - linpred;
    const double log_hazard = log_k - linpred + (k - 1.0) * log_scaled_time;
    const double log_surv = -std::exp(k * log_scaled_time);
    return log_hazard + log_surv;
  };
  auto log_weibull_tail_probability = [](double time, double linpred, double k) {
    const double log_time = std::log(time);
    const double log_scaled_time = log_time - linpred;
    return -std::exp(k * log_scaled_time);
  };
  auto dlog_weibull_density_dlinpred = [](double time, double linpred, double k) {
    const double log_time = std::log(time);
    const double log_scaled_time = log_time - linpred;
    const double exp_term = std::exp(k * log_scaled_time);
    return k * (exp_term - 1.0);
  };
  auto dlog_weibull_tail_dlinpred = [](double time, double linpred, double k) {
    const double log_time = std::log(time);
    const double log_scaled_time = log_time - linpred;
    const double exp_term = std::exp(k * log_scaled_time);
    return k * exp_term;
  };
  auto dlog_weibull_density_dlogk = [](double time, double linpred, double k) {
    const double log_time = std::log(time);
    const double log_scaled_time = log_time - linpred;
    const double exp_term = std::exp(k * log_scaled_time);
    const double ku = k * log_scaled_time;
    return 1.0 + ku - ku * exp_term;
  };
  auto dlog_weibull_tail_dlogk = [](double time, double linpred, double k) {
    const double log_time = std::log(time);
    const double log_scaled_time = log_time - linpred;
    const double exp_term = std::exp(k * log_scaled_time);
    const double ku = k * log_scaled_time;
    return -ku * exp_term;
  };

  const int n_event = event_time.n_elem;
  const int n_censor = censor_time.n_elem;
  const int n_cols = X_w_event.n_cols > 0 ? X_w_event.n_cols : X_w_censor.n_cols;
  if (model_type != 0 && model_type != 1) {
    Rcpp::stop("model_type must be 0 (exponential) or 1 (weibull).");
  }
  const bool is_weibull = model_type == 1;
  const int expected_params = is_weibull ? (n_cols + 2) : (n_cols + 1);
  if (params.n_elem != static_cast<arma::uword>(expected_params)) {
    Rcpp::stop("params length does not match model_type.");
  }
  if ((X_w_event.n_cols != static_cast<arma::uword>(n_cols)) ||
      (X_w_censor.n_cols != static_cast<arma::uword>(n_cols))) {
    Rcpp::stop("Event/censor design matrices must have the same number of columns.");
  }
  if (avg_linpred_event.n_elem != static_cast<arma::uword>(n_event) ||
      X_w_event.n_rows != static_cast<arma::uword>(n_event)) {
    Rcpp::stop("Event input lengths are inconsistent.");
  }
  if (avg_linpred_censor.n_elem != static_cast<arma::uword>(n_censor) ||
      X_w_censor.n_rows != static_cast<arma::uword>(n_censor)) {
    Rcpp::stop("Censor input lengths are inconsistent.");
  }

  const double beta_g = params[0];
  const arma::vec beta_w = params.subvec(1, n_cols);
  const double log_k = is_weibull ? params[n_cols + 1] : 0.0;
  const double k = is_weibull ? std::exp(log_k) : 1.0;

  const arma::vec xb_event = X_w_event * beta_w;
  const arma::vec xb_censor = X_w_censor * beta_w;

  double node_adj[kNumNodes];
  const double term_c = beta_g * post_c;
  for (int j = 0; j < kNumNodes; ++j) {
    node_adj[j] = term_c * kNodes[j];
  }

  const double* time_event_ptr = event_time.memptr();
  const double* avg_event_ptr = avg_linpred_event.memptr();
  const double* xb_event_ptr = xb_event.memptr();
  const double* x_event_ptr = X_w_event.memptr();

  const double* time_censor_ptr = censor_time.memptr();
  const double* avg_censor_ptr = avg_linpred_censor.memptr();
  const double* xb_censor_ptr = xb_censor.memptr();
  const double* x_censor_ptr = X_w_censor.memptr();

  arma::vec grad(params.n_elem, arma::fill::zeros);
  double total_ll = 0.0;

  auto accumulate_obs = [&](double time,
                            double avg_linpred,
                            double xb,
                            const double* x_ptr,
                            int n_rows,
                            int row_idx,
                            bool is_event,
                            arma::vec& grad_ref,
                            double& total_ll_ref) {
    double val[kNumNodes];
    double dlinpred[kNumNodes];
    double dlogk[kNumNodes];

    const double base = beta_g * avg_linpred + xb;
    double max_val = -1.0e300;

    for (int j = 0; j < kNumNodes; ++j) {
      const double linpred = base + node_adj[j];
      double log_f = 0.0;
      double dlogf_dlinpred = 0.0;
      double dlogf_dlogk = 0.0;

      if (!is_weibull) {
        if (is_event) {
          log_f = log_exp_density(time, linpred);
          dlogf_dlinpred = dlog_exp_density_dlinpred(time, linpred);
        } else {
          log_f = log_exp_tail_probability(time, linpred);
          dlogf_dlinpred = dlog_exp_tail_dlinpred(time, linpred);
        }
      } else {
        if (is_event) {
          log_f = log_weibull_density(time, linpred, log_k, k);
          dlogf_dlinpred = dlog_weibull_density_dlinpred(time, linpred, k);
          dlogf_dlogk = dlog_weibull_density_dlogk(time, linpred, k);
        } else {
          log_f = log_weibull_tail_probability(time, linpred, k);
          dlogf_dlinpred = dlog_weibull_tail_dlinpred(time, linpred, k);
          dlogf_dlogk = dlog_weibull_tail_dlogk(time, linpred, k);
        }
      }

      val[j] = log_f + kLogWeights[j];
      dlinpred[j] = dlogf_dlinpred;
      dlogk[j] = dlogf_dlogk;
      if (val[j] > max_val) {
        max_val = val[j];
      }
    }

    double sum_exp = 0.0;
    for (int j = 0; j < kNumNodes; ++j) {
      sum_exp += std::exp(val[j] - max_val);
    }

    total_ll_ref += max_val + std::log(sum_exp);
    const double inv_sum = 1.0 / sum_exp;
    
    double d_beta_w_sum = 0.0;

    for (int j = 0; j < kNumNodes; ++j) {
      const double weight = std::exp(val[j] - max_val) * inv_sum;
      const double dlinpred_weight = weight * dlinpred[j];
      const double dlinpred_dbetag = avg_linpred + post_c * kNodes[j];

      grad_ref[0] += dlinpred_weight * dlinpred_dbetag;
      
      if (is_weibull) {
        grad_ref[n_cols + 1] += weight * dlogk[j];
      }
      d_beta_w_sum += dlinpred_weight;
    }

    for (int col = 0; col < n_cols; ++col) {
      grad_ref[1 + col] += d_beta_w_sum * x_ptr[row_idx + n_rows * col];
    }
  };

  bool run_parallel = false;
#ifdef _OPENMP
  run_parallel = use_openmp;
#endif

  if (run_parallel) {
#pragma omp parallel
    {
      arma::vec grad_local(params.n_elem, arma::fill::zeros);
      double total_ll_local = 0.0;

#pragma omp for nowait
      for (int i = 0; i < n_event; ++i) {
        accumulate_obs(time_event_ptr[i],
                       avg_event_ptr[i],
                       xb_event_ptr[i],
                       x_event_ptr,
                       n_event,
                       i,
                       true,
                       grad_local,
                       total_ll_local);
      }

#pragma omp for nowait
      for (int i = 0; i < n_censor; ++i) {
        accumulate_obs(time_censor_ptr[i],
                       avg_censor_ptr[i],
                       xb_censor_ptr[i],
                       x_censor_ptr,
                       n_censor,
                       i,
                       false,
                       grad_local,
                       total_ll_local);
      }

#pragma omp critical
      {
        total_ll += total_ll_local;
        grad += grad_local;
      }
    }
  } else {
    for (int i = 0; i < n_event; ++i) {
      accumulate_obs(time_event_ptr[i],
                     avg_event_ptr[i],
                     xb_event_ptr[i],
                     x_event_ptr,
                     n_event,
                     i,
                     true,
                     grad,
                     total_ll);
    }

    for (int i = 0; i < n_censor; ++i) {
      accumulate_obs(time_censor_ptr[i],
                     avg_censor_ptr[i],
                     xb_censor_ptr[i],
                     x_censor_ptr,
                     n_censor,
                     i,
                     false,
                     grad,
                     total_ll);
    }
  }

  if (!R_finite(total_ll)) {
    return Rcpp::List::create(
      Rcpp::Named("value") = 1.0e12,
      Rcpp::Named("gradient") = arma::vec(params.n_elem, arma::fill::zeros)
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = -total_ll,
    Rcpp::Named("gradient") = -grad
  );
}