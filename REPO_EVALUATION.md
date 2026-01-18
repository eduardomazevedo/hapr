# Repository Evaluation

This review focuses on correctness, consistency, maintainability, and performance tradeoffs given the recent changes (MLE linear, MLE survival exponential/Weibull, removal of Cox).

## What looks solid

- **Clear two‑stage architecture**: `hapr_first_stage()` + `hapr_second_stage()` is clean and well encapsulated.
- **MLE path is separated and performance‑minded**: Rcpp likelihoods with precomputed constants, and data splitting for censoring is a good speed/clarity tradeoff.
- **Model type boundaries are explicit**: Two‑stage models remain `lm`/`probit`, MLE survival is separate with `model_type` switch for exponential/Weibull.
- **Testing infrastructure**: point‑estimate tests + artifacts + slow‑test gating are consistent and helpful.

## Problems / Risks

### 1) Roxygen / NAMESPACE drift risk
- **Issue**: Edits to docs and NAMESPACE depend on running `roxygen2::roxygenize()` to register S3 methods and `useDynLib`. If it’s forgotten, errors like missing `.Call` symbols can happen.
- **Impact**: Runtime failures after package rebuilds or installs.
- **Suggestion**: Add a CI check that compares NAMESPACE to the roxygen output or run roxygen in CI before tests.

### 2) MLE “delta” handling is inconsistent
- **Issue**: In MLE survival, `delta` is from optimizer and in LM `delta` is the log‑sigma. But prints and result structure still only loosely standardize this (`parameters$delta` is a vector but model meaning differs).
- **Impact**: Harder to program against MLE results in downstream code or compare across models.
- **Suggestion**: Document a standard convention for `parameters$delta` per model in code comments and docs (e.g., `log_sigma` for LM, `log_k` for Weibull). Optionally add a `delta_names` field or `delta_meaning` in result.

### 3) Stage‑1 and MLE print code depends on vcov naming
- **Issue**: MLE SE/CI printing relies on vcov dimnames. This was fixed for MLE by setting dimnames, but the pattern is fragile if optimizer returns unnamed params or users change `start_beta` names.
- **Impact**: Missing or mis‑aligned CI tables in prints.
- **Suggestion**: Enforce param naming prior to optimization and add a small internal validator to assert names line up (or set names on `opt$par` explicitly before Hessian inversion).

### 4) Potential mismatch in survival Weibull formulation
- **Issue**: Weibull model uses `scale = exp(linpred)` and shape `k = exp(log_k)`. The chosen log‑pdf/log‑survival formulas should be validated against the stated parameterization in `dev/theory.tex`.
- **Impact**: If the parameterization is off, estimates could be biased.
- **Suggestion**: Add a short comment block in C++ describing the exact PDF/S(t) parameterization; add a small unit test that compares loglik to base R `dweibull/pweibull` under a fixed parameterization.

### 5) Test helpers location (resolved)
- **Status**: Mock dataset helpers now live in `tests/testthat/helper-mock_dataset.R`, which is auto-loaded by testthat.

### 6) `hapr_simulate()` uses `fit$coefficients$theta`
- **Issue**: `hapr_simulate()` references `fit$coefficients$theta`, but `hapr_fit` stores parameters under `fit$parameters` not `fit$coefficients`.
- **Impact**: Potential runtime error or incorrect behavior depending on object structure.
- **Suggestion**: Update to `fit$parameters$theta` (or verify and standardize naming across outputs).

### 7) Print methods and summaries are growing in complexity
- **Issue**: Print logic is now somewhat dense and has many conditional paths.
- **Impact**: Harder to maintain and more likely to regress formatting.
- **Suggestion**: Consolidate printing into a small shared helper module (e.g., `R/print_helpers.R`) and keep per‑class print methods minimal.

### 8) First‑stage model type validation
- **Issue**: `hapr_first_stage()` now accepts `mle`, but several docs/man pages (after regeneration) might still mention only `lm`/`probit`.
- **Impact**: Confusion for users and inconsistent API contract.
- **Suggestion**: Ensure roxygen in `R/hapr_stage_1.R` fully documents `mle` and regenerate docs.

## Performance Opportunities

- **Vectorization vs branching**: The split‑event/censor strategy is good. Additional micro‑optimization likely has diminishing returns compared to BLAS tuning or fewer GH nodes.
- **Configurable quadrature nodes**: Allowing lower GH node counts could reduce runtime for quick iterations (e.g., a `quadrature_nodes` argument for MLE survival).
- **Potential caching**: Caching `X_w %*% beta_w` or node adjustments per iteration is already done in C++.

## Test Coverage Gaps

- **Weibull correctness vs known formula**: No direct test against analytical density/survival from base R.
- **MLE delta SE/CI**: Tests currently validate beta estimates, but not delta estimates or their SEs.
- **Edge cases**: No tests for all‑censored or all‑event cases; censoring split code should be validated in these extremes.

## General Suggestions

- **Document MLE limitations**: Explicitly state in docs that MLE SEs ignore first‑stage uncertainty.
- **Keep Rcpp exports stable**: Consider a short `tools/` script to run `compileAttributes()` and `roxygenize()` together before release.
- **Consistent object schema**: Ensure all fits (`hapr_fit`, `hapr_mle_fit`) store core parameters (`beta`, `theta`, `delta`) in consistent locations with standard names.

---

If you want, I can prioritize these into a short action plan or open issues by severity.
