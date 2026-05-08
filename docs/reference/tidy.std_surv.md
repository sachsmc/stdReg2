# Provide tidy output from a std_surv object for use in downstream computations

Tidy summarizes information about the components of the standardized
model fit.

## Usage

``` r
# S3 method for class 'std_surv'
tidy(x, ...)
```

## Arguments

- x:

  An object of class std_surv

- ...:

  Not currently used

## Value

A data.frame

## Examples

``` r
require(survival)
set.seed(8)
n <- 300
Z <- rnorm(n)
X <- rnorm(n, mean = Z)
time <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
U <- pmin(time, C) # time at risk
D <- as.numeric(time < C) # event indicator
dd <- data.frame(Z, X, U, D)
x <- standardize_coxph(
  formula = Surv(U, D) ~ X + Z + X * Z,
  data = dd, values = list(X = seq(-1, 1, 0.5)), times = c(2,3,4)
)

tidy(x)
#>         X   Estimate  Std.Error lower.0.95 upper.0.95 time contrast transform
#> X    -1.0 0.37425848 0.05636602 0.26378312  0.4847338    2     none  identity
#> X.1  -0.5 0.18100390 0.04207569 0.09853707  0.2634707    2     none  identity
#> X.2   0.0 0.11371977 0.02911956 0.05664648  0.1707931    2     none  identity
#> X.3   0.5 0.10190382 0.02572408 0.05148554  0.1523221    2     none  identity
#> X.4   1.0 0.10122665 0.02449415 0.05321900  0.1492343    2     none  identity
#> X1   -1.0 0.29008665 0.05502003 0.18224936  0.3979239    3     none  identity
#> X.11 -0.5 0.11993079 0.03651545 0.04836184  0.1914998    3     none  identity
#> X.21  0.0 0.08120193 0.02656267 0.02914006  0.1332638    3     none  identity
#> X.31  0.5 0.07969740 0.02411627 0.03243039  0.1269644    3     none  identity
#> X.41  1.0 0.08403060 0.02343156 0.03810559  0.1299556    3     none  identity
#> X2   -1.0 0.29008665 0.05502003 0.18224936  0.3979239    4     none  identity
#> X.12 -0.5 0.11993079 0.03651545 0.04836184  0.1914998    4     none  identity
#> X.22  0.0 0.08120193 0.02656267 0.02914006  0.1332638    4     none  identity
#> X.32  0.5 0.07969740 0.02411627 0.03243039  0.1269644    4     none  identity
#> X.42  1.0 0.08403060 0.02343156 0.03810559  0.1299556    4     none  identity
#>       measure
#> X    survival
#> X.1  survival
#> X.2  survival
#> X.3  survival
#> X.4  survival
#> X1   survival
#> X.11 survival
#> X.21 survival
#> X.31 survival
#> X.41 survival
#> X2   survival
#> X.12 survival
#> X.22 survival
#> X.32 survival
#> X.42 survival
```
