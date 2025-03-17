Improvements
Store objects more organized


Functions to implement
lm
    - hapr_lm(y, gc, w, improvement_ratio)
      - (w first stage, second stage options). Returns fit object.
    - simulate(fit, W, optional Gc)
      - Returns data frame with (W, Gc, Gf, Y). Same function for all models, except for Y.
      - Under the hood will need draw_gf function, common to all models.
    - predict(fit, W, Gc, Gf). G_c and G_f are optional. Returns Y_hat.



probit
  - Same thing. 
  - predict should have an option.




cox
  - same thing
  - Predict will give either linear predictor or risk ratio.
  - basehaz --> just mutiply the internal basehaz
  - survfit returns a survival curve given data.