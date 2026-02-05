# Plan: Reparameterize Survival MLE in Terms of Gamma

## Goal
Move the parametric survival MLE from estimating `beta` directly to estimating `gamma`, then map `gamma -> beta` with a delta-method SE. This yields correct SEs by accounting for first-stage uncertainty and simplifies the C++ likelihood.

## Key Idea (from `dev/theory.tex`)
- The model assumes hazard
  \[\lambda(t | \ell; \delta) = \lambda_0(t | \delta) \exp(\ell)\] where \(\ell = \beta_g G_f + \beta_w^\top W\).
- Conditional on observed \(G_c, W\), we can write
  \[G_f = a G_c + b W^\top \theta + c Z\] with \(Z \sim N(0,1)\).
- Lemma (expected hazard):
  \[\mathbb{E}[\lambda(t | L; \delta) | G_c, W] = e^{(\beta_g c)^2/2} \lambda_0(t|\delta) \exp(\gamma_g G_c + \gamma_w^\top W)\]
  where \(\gamma\) is the coefficient vector of a regression of \(Y\) on \(G_c, W\) (same linear mapping as in the linear model):
  \[\gamma_g = \beta_g a\]
  \[\gamma_w = \beta_g b \theta + \beta_w\].
- Hence the likelihood can be written in terms of \(\gamma\) (plus an adjustment factor involving \(\beta_g\)) and the hazard becomes a standard parametric survival likelihood with covariates \(G_c, W\).

## Mathematical Steps to Implement
1. **Define parameterization**
   - Replace \(\beta\) in the parametric survival likelihood with \(\gamma\).
   - Keep baseline hazard parameters \(\delta\) unchanged.

2. **Gamma-to-beta mapping** (linear conversion, same as linear model):
   - \(\beta_g = \gamma_g / a\)
   - \(\beta_w = \gamma_w - \beta_g b \theta\)

3. **Adjustment factor** in expected hazard:
   - Expected hazard contains multiplier \(\exp((\beta_g c)^2/2)\).
   - When \(\beta_g\) is expressed via \(\gamma_g\), this is
     \[\exp\left(\frac{c^2}{2} \left(\frac{\gamma_g}{a}\right)^2\right)\].
   - This term affects the baseline hazard scale.

4. **Likelihood with censoring**
   - For each observation, likelihood contribution is
     - Uncensored: \(f(t | \gamma, \delta)\)
     - Left-censored: \(F(t | \gamma, \delta)\)
     - Right-censored: \(1 - F(t | \gamma, \delta)\)
   - Use standard parametric survival density/CDF for exponential and Weibull with covariates \(\eta = \gamma_g G_c + \gamma_w^\top W\).

5. **Derive analytic gradients (C++)**
   - Likelihood wrt \(\gamma\) and \(\delta\) (baseline parameters).
   - Optional: keep existing gradient structure but substitute \(\eta\) and hazard multiplier.

6. **Delta-method standard errors for beta**
   - Let parameter vector estimated by MLE be \(\psi = (\gamma, \delta)\).
   - Compute Jacobian \(J = \partial \beta / \partial \gamma\).
   - Only \(\gamma\) affects \(\beta\); \(\delta\) block is zero.
   - For \(\beta_g = \gamma_g/a\): \(\partial \beta_g / \partial \gamma_g = 1/a\).
   - For \(\beta_w = \gamma_w - (b\theta)(\gamma_g/a)\):
     - \(\partial \beta_w / \partial \gamma_w = I\)
     - \(\partial \beta_w / \partial \gamma_g = -(b\theta)/a\)
   - Var(\(\beta\)) = \(J \, \text{Var}(\gamma) \, J^\top\).

## Implementation Plan
1. **Add new gamma-based MLE** (parallel implementation)
   - New R entry point and new C++ likelihood/gradient:
     - `R/hapr_mle_survival_gamma.R`
     - `src/hapr_mle_survival_gamma.cpp`
     - `src/hapr_mle_survival_gamma_grad.cpp`
   - Keep the old beta-based MLE unchanged.

2. **Wire into R**
   - Mirror function signature of existing survival MLE.
   - Store output with `parameters$gamma` and `parameters$delta`.
   - Provide `beta` in output by mapping gamma to beta (for parity).

3. **Tests for equivalence and runtime**
   - Create new tests comparing:
     - Point estimates from gamma-MLE vs beta-MLE (should match closely).
     - Log-likelihood values (should match).
     - Gradient norms at optimum (should be near zero and similar).
     - Runtime (gamma vs beta) using `system.time()` on identical data and seed.
   - Use simulated data (small n) with fixed seed.
   - Save artifacts to `tests/testthat/_artifacts` via `testthat::test_path()`:
     - Estimates and SEs for both methods.
     - Timing summaries.

4. **SE validation (expected direction)**
   - Verify that beta SEs from gamma-based MLE (delta-method) are larger than or equal to beta SEs from the old beta-based MLE, because the old method ignores theta uncertainty.
   - Record differences in artifacts for manual review.

5. **Migration**
   - After passing all tests, replace old MLE path with gamma-based path.
   - Keep old code behind a feature flag or remove after stabilization.

## Coding Notes
- Keep strict parameter ordering to align Jacobian with VCOV blocks.
- Maintain naming convention: `gc` and `gf`.
- Ensure all likelihood code uses normalized `G_c`.
- Prefer analytical gradients in C++ (matches current style).

## Next Steps
1. Locate existing survival MLE R/C++ code and mirror interfaces.
2. Implement gamma-based likelihood in C++.
3. Add R wrappers and delta-method SEs.
4. Write equivalence and runtime tests vs current implementation.
