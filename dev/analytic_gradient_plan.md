# Analytic gradient plan for survival MLE (exponential + Weibull)

## Goal
Add analytic gradients for the survival MLE log-likelihood (Gauss-Hermite quadrature) and integrate them into a new estimator while keeping the existing estimator intact. Validate against current likelihood and estimator with fast, artifact-producing tests and benchmarking.

## Scope
- Models: exponential and Weibull survival.
- Data: event time + censoring indicator (0/1), current PRS `gc`, covariates `w`.
- Likelihood: Gauss-Hermite integration already used in `hapr_mle_survival_nll_split_cpp`.
- Parameters: `beta_g`, `beta_w` (including intercept), and for Weibull `log_k` (shape).

## Mathematical derivations

### Notation
For each observation i, define
- `avg_i = a * gc_i + b * w_i' * theta`, from posterior mean of `G_f`.
- `xb_i = w_i' * beta_w` (including intercept in `w_i`).
- For node j: `z_j` with log weight `log_w_j`.
- `linpred_{ij} = beta_g * avg_i + xb_i + beta_g * post_c * z_j`.
- Log-sum-exp for event/censor contributions:
  - `log L_i = log( sum_j exp( log f_ij + log_w_j ) )`.

For gradients we use:
- `w_{ij} = exp( log f_ij + log_w_j - log L_i )` (softmax weights).
- `d log L_i / dθ = sum_j w_{ij} * d log f_ij / dθ`.

### Exponential model
Let `rate = exp(linpred)`.
- Event log density: `log f = linpred - rate * t`.
  - `d log f / d linpred = 1 - rate * t`.
- Censor log tail: `log f = -rate * t`.
  - `d log f / d linpred = -rate * t`.

Derivatives of `linpred`:
- `d linpred / d beta_g = avg_i + post_c * z_j`.
- `d linpred / d beta_w = w_i` (row vector).

### Weibull model
Let `k = exp(log_k)` and `u = log(t) - linpred`.
- Event log density:
  - `log f = log_k - linpred + (k - 1) * u - exp(k * u)`.
  - `d log f / d linpred = -k + k * exp(k * u)`.
  - `d log f / d log_k = 1 + k * u - k * u * exp(k * u)`.
- Censor log tail:
  - `log f = -exp(k * u)`.
  - `d log f / d linpred = k * exp(k * u)`.
  - `d log f / d log_k = -k * u * exp(k * u)`.

Chain rule for `beta_g` and `beta_w` is identical to the exponential case using `d log f / d linpred`.

## Implementation plan

### 1. New C++ likelihood with gradient
- Add `src/hapr_mle_survival_likelihood_grad.cpp` implementing:
  - `hapr_mle_survival_nll_split_grad_cpp(...)` that returns a list with `value` and `gradient`.
- Logic mirrors `hapr_mle_survival_nll_split_cpp`, but for each observation:
  - Compute per-node `log f_ij`.
  - Accumulate log-sum-exp for stability.
  - Compute softmax weights `w_{ij}`.
  - Accumulate gradients using formulas above.
- Return negative log-likelihood and **negative** gradient (since optimizer minimizes).

### 2. R wrapper for gradient-capable likelihood
- New R file (e.g. `R/hapr_mle_survival_likelihood_grad.R`).
- Provide `make_hapr_mle_likelihood_survival_grad(...)` that returns a list:
  - `fn(params)` returning scalar NLL.
  - `gr(params)` returning gradient vector.
- Use `hapr_mle_survival_nll_split_grad_cpp` under the hood.

### 3. New estimator using analytic gradient
- New estimator function (e.g. `hapr_mle_survival_grad`) analogous to `hapr_mle_survival`:
  - Same inputs, but uses `optim(..., gr = gr, method = "BFGS")`.
  - Keep outputs structured as `hapr_mle_fit` for consistency.
  - Keep original estimator unchanged.

## Testing + benchmarking plan

### 1. Likelihood/gradient unit tests (fast)
- New test file: `tests/testthat/test-hapr_mle_survival_grad_likelihood.R`.
- Generate small datasets (`n = 200`, two `var_epsilon` values, both model types).
- Compare:
  - Analytic gradient to finite-difference gradient (central diff) on a fixed parameter vector.
  - NLL values match old `hapr_mle_survival_nll_split_cpp` within tolerance.
- Produce artifact CSV with per-scenario:
  - Max abs grad diff, mean abs grad diff, NLL diff, timing for NLL + grad.

### 2. Estimator equivalence tests (fast)
- New test file: `tests/testthat/test-hapr_mle_survival_grad_estimator.R`.
- Compare `hapr_mle_survival` vs `hapr_mle_survival_grad` on small datasets:
  - Parameter estimates close (within tolerance).
  - Likelihood values close.
- Benchmark runtime of both estimators.
- Produce artifact CSV summarizing parameter diffs + runtimes.

### 3. Artifacts location
- Use `tests/testthat/_artifacts/mle_survival_grad/`.
- Write summary CSVs for easy inspection.

## Rollout steps
1. Implement new C++ likelihood + export in `RcppExports`.
2. Add R wrapper + new estimator function.
3. Add tests (likelihood/gradient + estimator equivalence).
4. Run targeted tests (survival-related) and verify artifacts.
5. Only proceed to additional refactors after results look good.
