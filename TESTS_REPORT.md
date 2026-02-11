# Test Inventory and Organization Review

## Test Inventory

- `tests/testthat/test-lowlevel_regression.R`
- `fit_lm` matches base `lm` with intercept (coefficients, vcov, sigma^2, R^2).
- `fit_lm` matches base `lm` without intercept.
- `fit_probit` matches `glm(..., family = binomial(link = "probit"))` with intercept (coefficients, vcov; R^2 bounded).
- `fit_probit` matches `glm` without intercept.
- `fit_lm` returns expected structure/types/dimensions and reasonable ranges.
- `fit_probit` returns expected structure/types/dimensions and reasonable ranges.
- Edge case: `fit_lm` with single predictor matches `lm`.
- Edge case: `fit_probit` with single predictor matches `glm`.

- `tests/testthat/test-hapr_stage_1_validation.R`
- `hapr_first_stage` rejects invalid `model_type` values.
- `hapr_first_stage` rejects `w` as data frame (expects numeric matrix).
- `hapr_first_stage` rejects constant columns in `w` (including all-zeros).
- `hapr_first_stage` rejects non–full-rank `w` (dependent columns).
- `hapr_first_stage` rejects empty `w` (0 columns).
- `hapr_first_stage` rejects mismatched dimensions across `y`, `gc`, `w`.
- `hapr_first_stage` rejects missing values in `y`, `gc`, or `w`.
- `hapr_first_stage` requires logical `y` for `model_type = "probit"` (rejects numeric 0/1).

- `tests/testthat/test-stage_2_point_estimates.R`
- For `model_type` in `{lm, probit}`, `n` in `{1e3, 1e4}` (plus `1e5` when `RUN_SLOW_TESTS=true`), and `var_epsilon` in `{0.5..0.9}`:
- Simulates data (`mock_dataset_lm`/`mock_dataset_probit`), runs `hapr_first_stage` + `hapr_second_stage`.
- Checks all coefficients (`gf`, intercept, `w1..w3`) are within 3 SE of truth.
- Writes per-scenario CSV artifacts and a summary CSV under `tests/testthat/_artifacts/point_estimates`.

- `tests/testthat/test-stage_2_coverage.R` (skipped unless `RUN_SLOW_TESTS=true`)
- For `model_type` in `{lm, probit}`, `n` in `{1e3, 1e4}`, `var_epsilon` in `{0.5..0.9}`:
- Runs 100 simulations per scenario; computes coverage and SE/SD ratio.
- Asserts coverage >= 85% for all coefficients and SE/SD ratio in [0.85, 2.0].
- Saves per-scenario summary CSVs and per-coefficient histograms under `tests/testthat/_artifacts/coverage` plus an overall summary CSV.
- Skips one probit scenario: `n=1e3`, `var_epsilon=0.9`.

- `tests/testthat/test-hapr_mle_survival_exp.R`
- For `model_type` in `{exponential, weibull}`, `n` in `{1e3, 1e4}` (plus `1e5` when `RUN_SLOW_TESTS=true`), `var_epsilon` in `{0.5..0.9}`, and `log_k` in `{-1, 0, 1}` (weibull only):
- Simulates data and runs `hapr_mle_survival`, checking all coefficients within 4 SE.
- Writes per-scenario CSV artifacts and a summary CSV under `tests/testthat/_artifacts/mle_survival`.
- Verifies `hapr_mle_survival` runs without OpenMP (`use_openmp = FALSE`).

- `tests/testthat/test-hapr_mle_survival_coverage.R` (skipped unless `RUN_SLOW_TESTS=true`)
- For `model_type` in `{exponential, weibull}`, `n` in `{1e3, 1e4}`, `var_epsilon` in `{0.5..0.9}`, and `log_k = 0` (weibull only):
- Runs 100 simulations per scenario; computes coverage and SE/SD ratio.
- Asserts coverage >= 85% for all coefficients and SE/SD ratio within bounds (upper bound loosened in one scenario).
- Skips/adjusts simulations where `improvement_ratio` exceeds `max_improvement_ratio` (tracked in meta CSV).
- Writes per-scenario summary and meta CSVs under `tests/testthat/_artifacts/coverage_survival` plus a summary CSV.

- `tests/testthat/test-zzz-aggregate-artifacts.R` (skipped unless `RUN_SLOW_TESTS=true`)
- Runs `dev/aggregate_artifacts.R` if present and asserts `coverage_summary.csv` is created under `tests/testthat/_artifacts`.

- `tests/testthat/helper-mock_dataset.R`
- Helpers for generating mock datasets used in stage-2 and survival tests (`lm`, `probit`, `exponential`, `weibull`).

## Organization Assessment

### What’s working
- Test coverage is broad: low-level regression correctness, input validation, two-stage point estimates, two-stage coverage, survival MLE point estimates, survival coverage, artifact aggregation.
- Slow tests are gated with `RUN_SLOW_TESTS=true` and produce useful diagnostics/artifacts.
- Coverage tests include detailed per-scenario logging and artifacts that aid debugging.

### Duplications / DRY concerns
- **Scenario loops are duplicated** across `test-stage_2_point_estimates.R` and `test-stage_2_coverage.R` with only small differences (sim count, assertions, artifact types).
- **Survival MLE loops** are duplicated between `test-hapr_mle_survival_exp.R` and `test-hapr_mle_survival_coverage.R` (scenario grid, data generation, naming, coefficient extraction, artifact layout).
- **Seed logic and scenario naming** are repeated in multiple tests, increasing the chance of accidental divergence.
- **Artifact directory creation** and summary CSV handling are repeated in several files with slightly different patterns.
- **Coefficient assembly / ordering** (`true_coef`, `estimates[names(true_coef)]`, etc.) repeats across test files.

### Opportunities to organize more elegantly

1. **Centralize scenario grids and iteration helpers**
- Create helper functions in `tests/testthat/helper-scenarios.R` (or extend `helper-mock_dataset.R`) that yield scenario grids for stage-2 and survival MLE tests.
- Example helpers:
  - `stage2_scenarios(run_slow)` -> list of scenario records (`model_type`, `n`, `var_epsilon`, `var_v`, `improvement_ratio`, `scenario_name`, `seed_base`).
  - `survival_scenarios(run_slow)` -> list with `model_type`, `n`, `var_epsilon`, `log_k`, `var_v`, `improvement_ratio`.
- This removes repeated loops and keeps grid definitions consistent.

2. **Shared functions for coefficient extraction / alignment**
- Helper functions that accept a `ci_beta` (or `fit`) and return aligned `estimates`, `se`, `ci` with a canonical order (`gf`, `(Intercept)`, `w1..w3`).
- Reduces repeated `names(true_coef)` alignment code and potential ordering mistakes.

3. **DRY artifact handling**
- A small artifact utility helper:
  - `ensure_artifact_dir(name)`
  - `write_scenario_csv(dir, scenario_name, df)`
  - `write_summary_csv(dir, name, df)`
- Helps avoid repeated `dir.exists`/`dir.create` boilerplate and keeps file naming uniform.

4. **Coverage vs point estimate tests as parameterized wrappers**
- For stage-2 tests, split into reusable functions:
  - `run_stage2_once(scenario)` -> returns fitted estimates and SEs.
  - `evaluate_point_estimates(scenario, result)` -> returns comparison table + boolean.
  - `evaluate_coverage(scenario, simulations)` -> returns summary table + coverage metrics.
- For survival tests, similarly split `run_survival_once` and evaluation functions.
- The test files then become short and declarative, with fewer repeated blocks.

5. **Consistent seeding strategy**
- Move seed calculation into helpers so seeds remain stable if scenario grids change.
- Make seed components explicit and reused between point and coverage tests to ease debugging.

6. **Unify logic for slow-test gating**
- Helper `is_slow_enabled()` (reads `RUN_SLOW_TESTS`) used across all tests to avoid repeated parsing.

### Suggested organization (minimal disruption)
- Add helper file(s):
  - `tests/testthat/helper-scenarios.R` (scenario grids + seed helpers)
  - `tests/testthat/helper-artifacts.R` (artifact dir + write helpers)
  - Optional: `tests/testthat/helper-estimates.R` (extract/alignment utils)
- Refactor stage-2 point + coverage tests to call shared helpers.
- Refactor survival point + coverage tests similarly.

This should significantly reduce duplication while keeping the tests readable and maintainable. If you want, I can propose specific helper function signatures and a refactor plan next.
