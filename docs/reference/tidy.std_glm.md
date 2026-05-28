# Provide tidy output from a std_glm object for use in downstream computations

Tidy summarizes information about the components of the standardized
regression fit.

## Usage

``` r
# S3 method for class 'std_glm'
tidy(x, ...)
```

## Arguments

- x:

  An object of class std_glm

- ...:

  Not currently used

## Value

A data.frame

## Examples

``` r
set.seed(6)
n <- 100
Z <- rnorm(n)
X <- rnorm(n, mean = Z)
Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
dd <- data.frame(Z, X, Y)
x <- standardize_glm(
  formula = Y ~ X * Z,
  family = "binomial",
  data = dd,
  values = list(X = 0:1),
  contrasts = c("difference", "ratio"),
  reference = 0
)
tidy(x)
#>   X   Estimate  Std.Error lower.0.95   upper.0.95   contrast transform
#> 1 0  0.5190639 0.06219793  0.3971582  0.640969573       none  identity
#> 2 1  0.3905311 0.09236927  0.2094907  0.571571540       none  identity
#> 3 0  0.0000000 0.00000000  0.0000000  0.000000000 difference  identity
#> 4 1 -0.1285328 0.06393487 -0.2538428 -0.003222729 difference  identity
#> 5 0  1.0000000 0.00000000  1.0000000  1.000000000      ratio  identity
#> 6 1  0.7523758 0.12876862  0.4999940  1.004757658      ratio  identity
```
