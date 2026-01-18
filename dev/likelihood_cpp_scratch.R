rm(list = ls())
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

# ==============================================================================
# 1. FINAL C++ IMPLEMENTATION
# ==============================================================================
sourceCpp(code = '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Hardcoded Gauss-Hermite Nodes (Pre-scaled: x * sqrt(2))
// These correspond to n=20 nodes for N(0,1) integration.
static const double NODES[20] = {
  -7.6190485, -6.5105902, -5.5787388, -4.7345813, -3.9439674, -3.1890148,
  -2.4586636, -1.7452473, -1.0429453, -0.3469642, 0.3469642, 1.0429453,
  1.7452473, 2.4586636, 3.1890148, 3.9439674, 4.7345813, 5.5787388,
  6.5105902, 7.6190485
};

// Hardcoded Log-Weights (Pre-scaled: log(w) - 0.5*log(pi))
static const double LOG_WEIGHTS[20] = {
  -29.7041538, -22.1163451, -16.6077588, -12.3335510, -8.9569762,
  -6.3033737, -4.2688523, -2.7886175, -1.8217647, -1.3440269,
  -1.3440269, -1.8217647, -2.7886175, -4.2688523, -6.3033737,
  -8.9569762, -12.3335510, -16.6077588, -22.1163451, -29.7041538
};

// [[Rcpp::export]]
double hapr_nll_final(const arma::vec& params,
                      const arma::vec& y,
                      const arma::vec& avg_linpred,   // Pre-mixed: (a*gc + b*w_theta)
                      const arma::mat& X_w,
                      double post_c) {
    
    // 1. Unpack Parameters
    double beta_g = params[0];
    int n_cols = X_w.n_cols;
    int n_obs = y.n_elem;
    
    // Subvectors (Zero-copy views)
    // Params structure: [beta_g, beta_w_1, ..., beta_w_k, delta]
    arma::vec beta_w = params.subvec(1, n_cols); 
    double delta = params[n_cols + 1];
    
    // 2. Pre-calculate Constants
    double sigma = std::exp(delta);
    double sigma_inv = 1.0 / sigma;
    double log_sigma = delta;
    
    // Matrix Multiplication (The only heavy linear algebra op)
    arma::vec xb = X_w * beta_w; 
    
    // Pre-calc node adjustments on stack (Scalar * Vector)
    double node_adj[20];
    double term_c = beta_g * post_c;
    for(int j=0; j<20; ++j) {
        node_adj[j] = term_c * NODES[j];
    }
    
    // 3. Raw Pointers (For maximum vectorization speed)
    const double* y_ptr = y.memptr();
    const double* avg_linpred_ptr = avg_linpred.memptr(); 
    const double* xb_ptr = xb.memptr();
    
    double total_ll = 0.0;
    
    // 4. Main Loop over Individuals
    for (int i = 0; i < n_obs; ++i) {
        
        // Base linear predictor: beta_g * (avg_linpred) + X*beta_w
        double base = beta_g * avg_linpred_ptr[i] + xb_ptr[i];
        
        double max_val = -1.0e20;
        double sum_exp = 0.0;
        
        // Inner Loop over Quadrature Nodes (Unrolled by compiler)
        for (int j = 0; j < 20; ++j) {
            double linpred = base + node_adj[j];
            double z = (y_ptr[i] - linpred) * sigma_inv;
            
            // Log PDF + Log Weight
            // log(dnorm) prop to -0.5 * z^2 - log(sigma)
            double val = (-0.5 * z * z - log_sigma) + LOG_WEIGHTS[j];
            
            // LogSumExp Logic
            if (val > max_val) {
                sum_exp = sum_exp * std::exp(max_val - val) + 1.0;
                max_val = val;
            } else {
                sum_exp += std::exp(val - max_val);
            }
        }
        total_ll += (max_val + std::log(sum_exp));
    }
    
    return -total_ll;
}
')

# ==============================================================================
# 2. R FACTORY
# ==============================================================================

#' Create HAPR Negative Log-Likelihood Function (Final C++ Version)
#'
#' @param y Outcome vector
#' @param gc Normalized PRS vector
#' @param w_theta Pre-calculated vector: w %*% theta
#' @param X_w Matrix of covariates including intercept
#' @param posterior List containing posterior constants (a, b, c)
#'
#' @return A function(params) compatible with optim
make_hapr_final <- function(y, gc, w_theta, X_w, posterior) {
  
  # OPTIMIZATION: Pre-mix the constant covariates
  # This vector is constant throughout the entire optimization.
  avg_linpred <- posterior$a * gc + posterior$b * w_theta
  
  # Return closure capturing the data
  function(params) {
    hapr_nll_final(params, y, avg_linpred, X_w, posterior$c)
  }
}

# ==============================================================================
# 3. VERIFICATION & TIMING
# ==============================================================================
set.seed(42)
N <- 10000 
cat(sprintf("Generating data for N = %d...\n", N))

y <- rnorm(N)
gc <- rnorm(N)
w_theta <- rnorm(N)
X_w <- cbind(1, matrix(rnorm(N * 3), N, 3)) 
posterior <- list(a = 0.8, b = 0.2, c = 0.5)

# Params: Beta_g (1), Beta_w (4), Delta (1)
test_params <- c(1.5, 0.5, 0.2, -0.1, 0.05, log(1.2)) 

# Create the function
fn_final <- make_hapr_final(y, gc, w_theta, X_w, posterior)

# Calculate Likelihood
ll_val <- fn_final(test_params)
cat(sprintf("Negative Log Likelihood: %.6f\n", ll_val))

cat("\nBenchmarking (1000 iterations)...\n")
# We expect this to be < 1ms per call
microbenchmark(
  Final = fn_final(test_params),
  times = 1000
)
