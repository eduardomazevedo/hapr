# HAPR: Heritability Adjusted Predictions for Polygenic Risk Scores
<!-- badges: start -->
<img src="man/figures/logo.svg" alt="HAPR hex logo" align="right" height="139" />
<!-- badges: end -->

There is an enormous amount of ongoing research on the practical utility of polygenic scores. But polygenic scores improve at a fast pace, so that these results become stale quickly. HAPR takes **currently available data** + **heritability estimates** to fit a statistical model of how a future polygenic score will perform.

Basic applications:
  - **Medicine**. A key question is to evaluate the clinical utility of polygenic risk scores, either on their own or in combination with other predictors. HAPR estimates predictive performance based both on currently available polygenic scores and on expected improvements as larger GWASs become available.
  - **Genomic social science**. Studies often estimate the effect of the true genetic predictor on a trait. Regressions using a currently available noisy polygenic risk score underestimate the true effect. HAPR estimates the true effect of the genetic predictor on the trait using current data combined with heritability estimates.
  - **Insurance economics**. Making genetic predictions available may change both pricing and consumer decisions in insurance markets. HAPR can be used to estimate the predictive power of both current polygenic risk scores and also how this power is expected to change as larger GWASs become available.

HAPR models quantitative traits (like height) and binary traits (like disease occurrence). HAPR uses robust, simple estimators that are scalable to large datasets.

## Installation

Install with

``` r
devtools::install_github("eduardomazevedo/hapr")
```

## Example
In this example, we'll analyze cardiovascular disease risk prediction. Our current GWAS polygenic score has an R^2 of 2%, whereas the heritability of the disease is about 20%. We can use HAPR to measure the effectiveness of the current polygenic score and a future score as the GWAS improves. The example considers a future score that reaches the R^2 of 20%.

``` r
library(hapr)
library(tidyverse)

# Sample data: electronic health records for cardiac patients
# - ascvd: Disease status (0/1)
# - polygenic_score: Current genetic risk score 
# - prevent_score: Traditional risk estimate
# - gender: Patient gender
head(heart_disease_data)
#> # A tibble: 6 × 4
#>   ascvd polygenic_score prevent_score gender
#>   <int>           <dbl>         <dbl> <fct> 
#> 1     0          -0.983        0.0428 female
#> 2     0          -0.797        0.181  male  
#> 3     1           0.749        0.0814 female
#> 4     0          -0.562        0.0147 male  
#> 5     0          -0.894        0.0238 male  
#> 6     0           0.577        0.0153 female

# Define covariates
covariates <- heart_disease_data |> select(prevent_score, gender)

# Run HAPR analysis
hapr_fit <- hapr(
  y = heart_disease_data$ascvd,          # Disease outcome
  gc = heart_disease_data$polygenic_score, # Current polygenic score
  w = covariates,                        # Covariates
  model_type = "probit",                 # Binary outcome model
  improvement_ratio = 10                 # Ratio of future R² (0.2) to current R² (0.02)
)

print(hapr_fit)
#> 🐵  HAPR (Heritability Adjusted Prediction) Model  🐵 
#> -------------------------------------------
#> Model type: probit 
#> 
#> Beta coefficients (future PRS effects):
#>               Estimate
#> (Intercept)    -2.5898
#> gf              2.6462
#> prevent_score   7.0468
#> gendermale      0.3171
#> 
#> Improvement ratio: 10.0000 ( r2_future )
#> R² current: 0.0200 ( user_provided )
#> R² future: 0.2000 
#> Max improvement ratio: 4573.9453
```

HAPR estimates the "true model" of disease risk based on both current data and estimated heritability. The data contains a current polygenic score `gc = polygenic_score`. The true model has a future polygenic score `gf`, which is what we would have if a larger GWAS reached the 20% heritability R^2. The estimated `beta` are the coefficients of a prediction model using `gf`. Under the hood, `HAPR` essentially assumes that the genes that will be discovered in the future correlate similarly to the covariates as the currently discovered genes. `HAPR` estimates the full model combining  feasible regressions with the current data and heritability estimates.

We can simulate what data would look like with the improved future polygenic score. This gives us a "dream dataset", as if we could know everyone's current and future prs, and everyone's estimated risk of disease based on models that use any combination of covariates, the current polygenic score, and future polygenic score.

``` r
# Simulate data with both current and future polygenic scores
simulated_data <- hapr_simulate(hapr_fit, covariates)

# Generate risk predictions
risk_dataset <- predict(hapr_fit, simulated_data)

# View the first few rows of predicted risk using either only covariates (w),
# covariates and the current polygenic score and covariates (gc_w)
# and risk using the future polygenic score and covariates (gf_w)
head(risk_dataset)
#>           gf         gc prevent_score gender    y_hat_w y_hat_gc_w  y_hat_gf_w
#> 1 -0.2303620 -1.0964639    0.04275714 female 0.04831298 0.02639364 0.001877092
#> 2 -0.1819389 -1.0103211    0.18114824   male 0.21443238 0.15604899 0.069752551
#> 3 -0.1400503  0.6446940    0.08136762 female 0.07629351 0.08999877 0.008492222
#> 4  0.1037776  0.6331684    0.01474663   male 0.03699994 0.04588312 0.029101187
#> 5  0.1831262  0.8419155    0.02379530   male 0.04159394 0.05617229 0.052569730
#> 6  0.5425069  1.8761863    0.01532521 female 0.03395948 0.06825223 0.147726142
```

See the vignettes for more examples with quantitative traits and more advanced analysis.


## Current repo structure

The package's core function is `hapr`.

Implementation uses a two-stage estimation approach:

1. First stage (`hapr_first_stage`): Runs feasible regressions using available data, without incorporating the improvement ratio (expected future R² / current R²)

2. Second stage (`hapr_second_stage`): Uses the improvement ratio to estimate the full model parameters

This structure allows easy calculation of results under different heritability assumptions by varying the improvement ratio parameter in the second stage, without having to redo the first stage, which has the bulk of computations.

## Current functionality
- Implemented models: lm, probit.
- Estimation (`hapr`), simulation (`hapr_simulate`), prediction (`predict`).

## TODOs:
- Implement standard errors with the delta method.
- Create coverage tests.
- Create more unit tests.
- Better vignettes / readme.
- General improvements.
- Fix a bunch of little issues required to submit to CRAN.
