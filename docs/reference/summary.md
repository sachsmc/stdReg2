# Summarizes parfrailty fit

This is a `summary` method for class `"parfrailty"`.

## Usage

``` r
# S3 method for class 'parfrailty'
summary(
  object,
  ci_type = "plain",
  ci_level = 0.95,
  digits = max(3L, getOption("digits") - 3L),
  ...
)
```

## Arguments

- object:

  an object of class `"parfrailty"`.

- ci_type:

  string, indicating the type of confidence intervals. Either "plain",
  which gives untransformed intervals, or "log", which gives
  log-transformed intervals.

- ci_level:

  desired coverage probability of confidence intervals, in decimal form.

- digits:

  the number of significant digits to use when printing..

- ...:

  not used.

## Value

An object of class "summary.parfrailty", which is a list that contains
relevant summary statistics about the fitted model

## See also

[`parfrailty`](https://sachsmc.github.io/stdReg2/reference/parfrailty.md)

## Author

Arvid Sjölander and Elisabeth Dahlqwist.

## Examples

``` r
## See documentation for parfrailty
```
