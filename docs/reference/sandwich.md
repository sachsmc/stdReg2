# Compute the sandwich variance components from a model fit

Compute the sandwich variance components from a model fit

## Usage

``` r
sandwich(fit, data, weights, t, fit.detail)
```

## Arguments

- fit:

  A fitted model object of class
  [glm](https://rdrr.io/r/stats/glm.html),
  [coxph](https://rdrr.io/pkg/survival/man/coxph.html), ah, or
  [survfit](https://rdrr.io/pkg/survival/man/survfit.html)

- data:

  The data used to fit the model

- weights:

  Optional weights

- t:

  Optional fixed time point for survival objects

- fit.detail:

  For Cox models, the result of running
  [coxph.detail](https://rdrr.io/pkg/survival/man/coxph.detail.html) on
  the model fit

## Value

A list consisting of the Fisher information matrix (I) and the Score
equations (U)
