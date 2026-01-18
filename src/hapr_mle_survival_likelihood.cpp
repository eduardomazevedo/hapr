#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
double hapr_mle_survival_exp_nll_cpp(const arma::vec& params,
                                     const arma::vec& event_time,
                                     const arma::vec& event_status,
                                     const arma::vec& avg_linpred,
                                     const arma::mat& X_w,
                                     double post_c) {
  auto log_exp_density = [](double time, double linpred) {
    const double rate = std::exp(linpred);
    return linpred - rate * time;
  };
  auto log_exp_tail_probability = [](double time, double linpred) {
    const double rate = std::exp(linpred);
    return -rate * time;
  };

  const int n_obs = event_time.n_elem;
  const int n_cols = X_w.n_cols;
  if (params.n_elem != static_cast<arma::uword>(n_cols + 1)) {
    Rcpp::stop("params length must equal ncol(X_w) + 1.");
  }
  if (event_status.n_elem != static_cast<arma::uword>(n_obs) ||
      avg_linpred.n_elem != static_cast<arma::uword>(n_obs) ||
      X_w.n_rows != static_cast<arma::uword>(n_obs)) {
    Rcpp::stop("Input lengths are inconsistent.");
  }

  const double beta_g = params[0];
  const arma::vec beta_w = params.subvec(1, n_cols);

  const arma::vec xb = X_w * beta_w;

  double node_adj[kNumNodes];
  const double term_c = beta_g * post_c;
  for (int j = 0; j < kNumNodes; ++j) {
    node_adj[j] = term_c * kNodes[j];
  }

  const double* time_ptr = event_time.memptr();
  const double* status_ptr = event_status.memptr();
  const double* avg_ptr = avg_linpred.memptr();
  const double* xb_ptr = xb.memptr();

  double total_ll = 0.0;
  for (int i = 0; i < n_obs; ++i) {
    const double base = beta_g * avg_ptr[i] + xb_ptr[i];

    double max_val = -1.0e300;
    double sum_exp = 0.0;

    for (int j = 0; j < kNumNodes; ++j) {
      const double linpred = base + node_adj[j];
      const double log_pdf = log_exp_density(time_ptr[i], linpred);
      const double log_surv = log_exp_tail_probability(time_ptr[i], linpred);
      const double val = status_ptr[i] * log_pdf +
        (1.0 - status_ptr[i]) * log_surv + kLogWeights[j];

      if (val > max_val) {
        sum_exp = sum_exp * std::exp(max_val - val) + 1.0;
        max_val = val;
      } else {
        sum_exp += std::exp(val - max_val);
      }
    }

    total_ll += max_val + std::log(sum_exp);
  }

  if (!R_finite(total_ll)) {
    return 1.0e12;
  }
  return -total_ll;
}

// [[Rcpp::export]]
double hapr_mle_survival_nll_split_cpp(const arma::vec& params,
                                       const arma::vec& event_time,
                                       const arma::vec& avg_linpred_event,
                                       const arma::mat& X_w_event,
                                       const arma::vec& censor_time,
                                       const arma::vec& avg_linpred_censor,
                                       const arma::mat& X_w_censor,
                                       double post_c,
                                       int model_type) {
  auto log_exp_density = [](double time, double linpred) {
    const double rate = std::exp(linpred);
    return linpred - rate * time;
  };
  auto log_exp_tail_probability = [](double time, double linpred) {
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

  const double* time_censor_ptr = censor_time.memptr();
  const double* avg_censor_ptr = avg_linpred_censor.memptr();
  const double* xb_censor_ptr = xb_censor.memptr();

  double total_ll = 0.0;

  if (!is_weibull) {
    for (int i = 0; i < n_event; ++i) {
      const double base = beta_g * avg_event_ptr[i] + xb_event_ptr[i];
      double max_val = -1.0e300;
      double sum_exp = 0.0;

      for (int j = 0; j < kNumNodes; ++j) {
        const double linpred = base + node_adj[j];
        const double val = log_exp_density(time_event_ptr[i], linpred) +
          kLogWeights[j];

        if (val > max_val) {
          sum_exp = sum_exp * std::exp(max_val - val) + 1.0;
          max_val = val;
        } else {
          sum_exp += std::exp(val - max_val);
        }
      }

      total_ll += max_val + std::log(sum_exp);
    }

    for (int i = 0; i < n_censor; ++i) {
      const double base = beta_g * avg_censor_ptr[i] + xb_censor_ptr[i];
      double max_val = -1.0e300;
      double sum_exp = 0.0;

      for (int j = 0; j < kNumNodes; ++j) {
        const double linpred = base + node_adj[j];
        const double val = log_exp_tail_probability(time_censor_ptr[i], linpred) +
          kLogWeights[j];

        if (val > max_val) {
          sum_exp = sum_exp * std::exp(max_val - val) + 1.0;
          max_val = val;
        } else {
          sum_exp += std::exp(val - max_val);
        }
      }

      total_ll += max_val + std::log(sum_exp);
    }
  } else {
    for (int i = 0; i < n_event; ++i) {
      const double base = beta_g * avg_event_ptr[i] + xb_event_ptr[i];
      double max_val = -1.0e300;
      double sum_exp = 0.0;

      for (int j = 0; j < kNumNodes; ++j) {
        const double linpred = base + node_adj[j];
        const double val = log_weibull_density(time_event_ptr[i], linpred, log_k, k) +
          kLogWeights[j];

        if (val > max_val) {
          sum_exp = sum_exp * std::exp(max_val - val) + 1.0;
          max_val = val;
        } else {
          sum_exp += std::exp(val - max_val);
        }
      }

      total_ll += max_val + std::log(sum_exp);
    }

    for (int i = 0; i < n_censor; ++i) {
      const double base = beta_g * avg_censor_ptr[i] + xb_censor_ptr[i];
      double max_val = -1.0e300;
      double sum_exp = 0.0;

      for (int j = 0; j < kNumNodes; ++j) {
        const double linpred = base + node_adj[j];
        const double val = log_weibull_tail_probability(time_censor_ptr[i], linpred, k) +
          kLogWeights[j];

        if (val > max_val) {
          sum_exp = sum_exp * std::exp(max_val - val) + 1.0;
          max_val = val;
        } else {
          sum_exp += std::exp(val - max_val);
        }
      }

      total_ll += max_val + std::log(sum_exp);
    }
  }

  if (!R_finite(total_ll)) {
    return 1.0e12;
  }
  return -total_ll;
}
