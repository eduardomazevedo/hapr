# HAPR: Heritability Adjusted Predictions for Polygenic Risk Scores
<!-- badges: start -->
<img src="man/figures/logo.svg" alt="HAPR hex logo" align="right" height="139" />
<!-- badges: end -->

There is an enormous amount of ongoing research on the practical utility of polygenic scores. But polygenic scores improve at a fast pace, so that these results become stale quickly. HAPR takes **currently available data** + **heritability estimates** to fit a statistical model of how a future polygenic score will perform.

Basic applications:

  - **Medicine**. A key question is to evaluate the clinical utility of polygenic risk scores, either on their own or in combination with other predictors. HAPR estimates predictive performance based both on currently available polygenic scores and on expected improvements as larger GWASs become available.
  
  - **Genomic social science**. Studies often estimate the effect of the true genetic predictor on a trait. Regressions using a currently available noisy polygenic risk score underestimate the true effect. HAPR estimates the true effect of the genetic predictor on the trait using current data combined with heritability estimates.

  - **Insurance economics**. Making genetic predictions available may change both pricing and consumer decisions in insurance markets. HAPR can be used to estimate the predictive power of both current polygenic risk scores and also how this power is expected to change as larger GWASs become available.

HAPR supports two-stage regression-based estimators for quantitative traits (linear) and binary traits (probit). It also provides two-stage MLE estimators for parametric survival models (exponential, Weibull).

## Installation

Install with

``` r
devtools::install_github("eduardomazevedo/hapr")
```

## Example
In this example, we'll analyze cardiovascular disease risk prediction. Our current GWAS polygenic score has an R^2 around 2-3%, whereas the heritability of the disease is about 25-30%. We can use HAPR to measure the effectiveness of the current polygenic score and a future score as the GWAS improves. The example considers a future score with a 10x improvement ratio.

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
# HAPR expects w as a numeric matrix without an intercept column.
w <- model.matrix(~ prevent_score + gender, data = covariates)[, -1, drop = FALSE]

# Run HAPR analysis
hapr_fit <- hapr(
  y = heart_disease_data$ascvd,          # Disease outcome
  gc = heart_disease_data$polygenic_score, # Current polygenic score
  w = w,                                 # Covariates (numeric matrix, no intercept)
  model_type = "probit",                 # Binary outcome model
  improvement_ratio = 10                 # Ratio of future R² to current R²
)

print(hapr_fit)
#> 🐵  HAPR (Heritability Adjusted Prediction) Model  🐵 
#> -------------------------------------------
#> Model type: probit 
#> 
#> Beta coefficients (future PRS effects):
#>             Estimate Std.Error   Lower   Upper
#> gf            2.6462    1.9631 -1.2013  6.4937
#> (Intercept)  -2.5898    0.7013 -3.9644 -1.2153
#> w1            7.0468    2.1548  2.8235 11.2702
#> w2            0.3171    0.3166 -0.3034  0.9376
#> 
#> Note: Standard errors are delta-method approximations and may be conservative.
#> 
#> Theta coefficients (gc ~ w): 
#>             Estimate Std.Error    Lower   Upper
#> (Intercept)  0.03132   0.04856 -0.06385 0.12649
#> w1           0.31124   0.53488 -0.73710 1.35958
#> w2          -0.09909   0.06672 -0.22987 0.03169
#> 
#> Stage 1 variance (v + epsilon): 
#>                        Estimate Std.Error Lower Upper
#> var_v_plus_var_epsilon   0.9998   0.04478 0.912 1.088
#> 
#> Improvement ratio: 10.0000 
#> R² current: 0.0274 
#> R² future: 0.2737 
#> Max improvement ratio: 4573.9453
```

HAPR estimates the "true model" of disease risk based on both current data and estimated heritability. The data contains a current polygenic score `gc = polygenic_score`. The true model has a future polygenic score `gf`, which is what we would have if a larger GWAS reached the 20% heritability R^2. The estimated `beta` are the coefficients of a prediction model using `gf`. Under the hood, `HAPR` essentially assumes that the genes that will be discovered in the future correlate similarly to the covariates as the currently discovered genes. `HAPR` estimates the full model combining  feasible regressions with the current data and heritability estimates.

We can simulate what data would look like with the improved future polygenic score. This gives us a "dream dataset", as if we could know everyone's current and future prs, and everyone's estimated risk of disease based on models that use any combination of covariates, the current polygenic score, and future polygenic score.

``` r
# Simulate data with both current and future polygenic scores
simulated_data <- hapr_simulate(hapr_fit, covariates)

# Build w in the same format used by hapr()
w_pred <- model.matrix(~ prevent_score + gender, data = simulated_data)[, -1, drop = FALSE]

risk_dataset <- predict(
  hapr_fit,
  list(w = w_pred, gc = simulated_data$gc, gf = simulated_data$gf)
)

head(risk_dataset)
#>   w.prevent_score w.gendermale         gc         gf    y_hat_w y_hat_gc_w
#> 1      0.04275714            0 -1.0964639 -0.2303620 0.04831298 0.02639364
#> 2      0.18114824            1 -1.0103211 -0.1819389 0.21443238 0.15604899
#> 3      0.08136762            0  0.6446940 -0.1400503 0.07629351 0.08999877
#> 4      0.01474663            1  0.6331684  0.1037776 0.03699994 0.04588312
#> 5      0.02379530            1  0.8419155  0.1831262 0.04159394 0.05617229
#> 6      0.01532521            0  1.8761863  0.5425069 0.03395948 0.06825223
#>    y_hat_gf_w
#> 1 0.001877092
#> 2 0.069752551
#> 3 0.008492222
#> 4 0.029101187
#> 5 0.052569730
#> 6 0.147726142
```

See the vignettes for more examples with quantitative traits and more advanced analysis.


## Functionality
Two-stage estimators can be run with `hapr()` or with `hapr_first_stage()` + `hapr_second_stage()`.
- Two-stage regression-based estimators: lm, probit.
- Two-stage MLE estimators: survival (exponential, Weibull).
- Estimation (`hapr`, `hapr_mle_survival`), simulation (`hapr_simulate`), prediction (`predict`).
