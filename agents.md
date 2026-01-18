# HAPR Development Guide

## Overview

HAPR (Heritability Adjusted Prediction) is an R package that estimates statistical models using future polygenic risk scores (PRS) that aren't yet available, by combining current PRS data with heritability estimates. The theoretical foundation is described in `dev/theory.tex`.

## Core Concept

- **G_c** (current PRS): Observed noisy polygenic risk score from current GWAS
- **G_f** (future PRS): True/unobserved genetic predictor that would be available from a larger future GWAS
- **Improvement ratio**: Ratio of R²_future / R²_current (expected improvement in PRS predictive power)
- **Beta coefficients**: True effects of G_f on outcome Y (what we want to estimate)

## Architecture: Two-Stage Estimation

The package uses a computationally efficient two-stage approach:

1. **First stage** (`hapr_first_stage`): Runs feasible regressions with current data:
   - `G_c ~ W` (theta coefficients)
   - `Y ~ G_c + W` (gamma coefficients)
   - These don't require the improvement ratio

2. **Second stage** (`hapr_second_stage`): Converts gamma to beta using the improvement ratio:
   - Calculates posterior parameters (a, b, c) from improvement ratio
   - Applies model-specific formulas to derive beta from gamma
   - Computes standard errors via delta method (analytical Jacobian)

This design allows exploring different improvement ratios without re-running the computationally expensive first stage.

## Supported Models

- **Two-stage estimators**: `lm` (linear) and `probit`
- **MLE estimators**: `lm` (development) and survival (`exponential`, `weibull`)

Model-specific conversion formulas are in `calculate_parameters()` and `calculate_analytical_jacobian()` for two-stage estimators.

## Key Files

- `R/hapr_wrapper.R`: Main `hapr()` function (combines stages 1 & 2)
- `R/hapr_stage_1.R`: First stage estimation
- `R/hapr_stage_2.R`: Second stage estimation, parameter conversion, delta method SEs
- `R/lowlevel_regression.R`: Core regression functions (`fit_lm`, `fit_probit`)
- `R/hapr_mle_lm.R`: MLE for linear model (development/testing)
- `R/hapr_mle_survival.R`: MLE for parametric survival models
- `src/hapr_mle_likelihood.cpp`: C++ linear MLE likelihood
- `src/hapr_mle_survival_likelihood.cpp`: C++ survival MLE likelihoods
- `R/preprocess.R`: Data preprocessing and validation
- `R/predict.R`: Prediction methods for `hapr_fit` objects
- `R/simulate.R`: Simulates data with both G_c and G_f
- `R/abc.R`: Helper to calculate posterior parameters a, b, c

## Development Notes

- **Working directory**: All scripts and commands are run from the project root directory
- **Assumption**: Always run commands from the project root; tests and scripts assume root paths
- **Normalization**: `G_c` must be normalized (unit variance) - checked in preprocessing
- **Coefficient naming**: Uses "gc" for current PRS, "gf" for future PRS in coefficient vectors
- **Delta method**: Standard errors computed via analytical Jacobian (see `calculate_analytical_jacobian`)
- **Testing**: Uses `testthat` - see `tests/testthat/` directory
- **Local test runs**: Prefer `devtools::test()` to ensure local sources are loaded
- **Test artifacts**: Write generated outputs to `tests/testthat/_artifacts` using `testthat::test_path()`
- **Slow tests**: Gate large `n` scenarios behind `RUN_SLOW_TESTS=true`

## Theory References

- Linear model: Section 5.1 in `dev/theory.tex`
- Probit model: Section 5.2 (formulas from @azevedo2024genetic)
- General theory: Sections 1-4 (model assumptions, identification, likelihood)

## Common Tasks

**Adding a new model type:**
1. Add model-specific conversion in `calculate_parameters()`
2. Add Jacobian calculation in `calculate_analytical_jacobian()`
3. Add low-level regression function in `lowlevel_regression.R`
4. Update `hapr_first_stage()` to handle new type

**Modifying standard errors:**
- Edit `calculate_analytical_jacobian()` to adjust derivatives
- Ensure column/row names match between Jacobian and covariance matrices

**Testing:**
- Unit tests in `tests/testthat/`
- Snapshots in `tests/testthat/_snaps/`
- Run with `devtools::test()` or `testthat::test_file()`
