
# HAPR: Heritability Adjusted Predictions for Polygenic Risk Scores

<!-- badges: start -->
<img src="assets/logo.svg" align="right" height="139" />
<!-- badges: end -->

Polygenic risk score predictions improve constantly, so that results become constantly out of date. HAPR estimates predictive performance based on expected GWAS improvements given heritability estimates.

Basic applications:
  - **Medicine**. A key question is to evaluate the clinical utility of polygenic risk scores, either on their own or in combination with other predictors. HAPR estimates predictive performance based both on currently available polygenic scores but also on expected improvements as larger GWASs become available.
  - **Genomic social science**. Studies often aim to calculate the effect of the true genetic predictor on a trait. Regressions using the currently availale noisy poligenic risk scores understimate the true effect. HAPR can be used to estimate the true effect of the genetic predictor on the trait.
  - **Insurance economics**. Making genetic predictions available may change both pricing and consumer decisions in insurance markets. HAPR can be used to estimate the predictive power of both current polygenic risk scores and also how this power is expected to change as larger GWASs become available.

HAPR models quantitative traits (like height), binary traits (like disease occurance), and survival data (like time of disease occurance). HAPR uses robust, simple estimators that are scalable to large datasets.

## Installation

You can install the development version of hapr like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example
Fit the model:

``` r
library(hapr)

height_fit <- hapr_lm(y = height, gc = height_prs, w = covariates, r2_current = 0.4, r2_future = 0.6)>

cad_fit <- hapr_probit(y = coronary_artery_disease, gc = cad_prs, w = covariates, r2_current = 0.03, r2_future = 0.20)>

cad_survival_fit <- hapr_cox(t = event_time, y = coronary_artery_disease, gc = cad_prs, w = covariates, r2_current = 0.03, r2_future = 0.20)>
```

TODO: make some cool plot

## Current repo structure

    - Three functions
    - hapr_lm, hapr_probit, hapr_cox.
    - Meat of the functions is in hapr_lm_first_stage, hapr_lm_second_stage, etc. We use the two stage estimators form the pdf. The first stage runs some regressions / standard statistical models. The second stage puts it together given heritability estimates.
    - The first stage runs the regressions. Takes as inputs y, gc, w. Does not need an improvement ratio or future r2. Returns a first-stage fit object.
    - The second stage takes the output of the first stage and an improvement ratio as parameters. And estimates the full model. Returns a fit object.
    - This structure makes it easy to calculate the results under different heritability estimates, since most of the computation happens in the first stage.
    - Helpter functions
      - regression functions to be run in the first stage.
      - preprocess_inputs.
      - abc to calculate abc constants.

## Current functionality
    - hapr_lm, hapr_probit work and I tested point estimates.
    - hapr_cox works, but I have not tested point estimates. The base hazard is not calculated yet.

## Future plans:
    - confidence intervals
    - coverage simulations
    - printing functions
    - docs
    - functions to create simulated datasets after fitting the model.
    - plotting functions