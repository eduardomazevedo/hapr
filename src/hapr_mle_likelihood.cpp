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
double hapr_mle_lm_nll_cpp(const arma::vec& params,
                           const arma::vec& y,
                           const arma::vec& avg_linpred,
                           const arma::mat& X_w,
                           double post_c) {
  const int n_obs = y.n_elem;
  const int n_cols = X_w.n_cols;
  if (params.n_elem != static_cast<arma::uword>(n_cols + 2)) {
    Rcpp::stop("params length must equal ncol(X_w) + 2.");
  }

  const double beta_g = params[0];
  const arma::vec beta_w = params.subvec(1, n_cols);
  const double log_sigma = params[n_cols + 1];
  const double sigma_inv = std::exp(-log_sigma);

  const arma::vec xb = X_w * beta_w;

  double node_adj[kNumNodes];
  const double term_c = beta_g * post_c;
  for (int j = 0; j < kNumNodes; ++j) {
    node_adj[j] = term_c * kNodes[j];
  }

  const double* y_ptr = y.memptr();
  const double* avg_ptr = avg_linpred.memptr();
  const double* xb_ptr = xb.memptr();

  double total_ll = 0.0;
  for (int i = 0; i < n_obs; ++i) {
    const double base = beta_g * avg_ptr[i] + xb_ptr[i];

    double max_val = -1.0e300;
    double sum_exp = 0.0;

    for (int j = 0; j < kNumNodes; ++j) {
      const double linpred = base + node_adj[j];
      const double z = (y_ptr[i] - linpred) * sigma_inv;
      const double val = (-0.5 * z * z - log_sigma) + kLogWeights[j];

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
