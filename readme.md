
# HAPR: Heritability Adjusted Predictions for Polygenic Risk Scores
<!-- badges: start -->
<img src="assets/logo.svg" align="right" height="139" />
<!-- badges: end -->

Polygenic risk score predictions improve constantly, so that results become constantly out of date. HAPR estimates predictive performance based on expected GWAS improvements given heritability estimates.

Basic applications:
  - **Medicine**. A key question is to evaluate the clinical utility of polygenic risk scores, either on their own or in combination with other predictors. HAPR estimates predictive performance based both on currently available polygenic scores but also on expected improvements as larger GWASs become available.
  - **Genomic social science**. Studies often estimate the effect of the true genetic predictor on a trait. Regressions using a currently availale noisy poligenic risk score understimate the true effect. HAPR estimates the true effect of the genetic predictor on the trait using current data combined with heritability estimates.
  - **Insurance economics**. Making genetic predictions available may change both pricing and consumer decisions in insurance markets. HAPR can be used to estimate the predictive power of both current polygenic risk scores and also how this power is expected to change as larger GWASs become available.

HAPR models quantitative traits (like height), binary traits (like disease occurance), and survival data (like time of disease occurance). HAPR uses robust, simple estimators that are scalable to large datasets.

## Installation

Install with

``` r
devtools::install_github("eduardomazevedo/hapr")
```

## Example
Fit the model:

``` r
library(hapr)

height_fit <- hapr_lm(y = height, gc = height_prs, w = covariates, r2_current = 0.4, r2_future = 0.6)>

cad_fit <- hapr_probit(y = coronary_artery_disease, gc = cad_prs, w = covariates, r2_current = 0.03, r2_future = 0.20)>

cad_survival_fit <- hapr_cox(t = event_time, y = coronary_artery_disease, gc = cad_prs, w = covariates, r2_current = 0.03, r2_future = 0.20)>
```

TODO: add a cool plot.

## Current repo structure

The package has three core functions:
- `hapr_lm`: For linear models
- `hapr_probit`: For binary outcomes 
- `hapr_cox`: For survival analysis

The implementation uses a two-stage estimation approach:

1. First stage (`hapr_*_first_stage`): Runs feasible regressions using available data, without incorporating the improvement ratio (expected future R² / current R²)

2. Second stage (`hapr_*_second_stage`): Uses the improvement ratio to estimate the full model parameters

This structure allows easy calculation of results under different heritability assumptions by varying the improvement ratio parameter.

## Current functionality
Core functions hapr_lm, hapr_probit, hapr_cox are working.

## Future plans:
   - output baseline hazard in hapr_cox.
   - confidence intervals.
   - coverage simulations.
   - printing functions.
   - functions to create simulated datasets after fitting the model.
   - plotting functions. 
   - docs.