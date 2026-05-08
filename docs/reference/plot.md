# Plots regression standardization fit

This is a `plot` method for class `"std_glm"`.

## Usage

``` r
# S3 method for class 'std_glm'
plot(
  x,
  plot_ci = TRUE,
  ci_type = "plain",
  ci_level = 0.95,
  transform = NULL,
  contrast = NULL,
  reference = NULL,
  summary_fun = "summary_std_glm",
  ...
)
```

## Arguments

- x:

  An object of class `"std_glm"`.

- plot_ci:

  if `TRUE`, add the confidence intervals to the plot.

- ci_type:

  A string, indicating the type of confidence intervals. Either "plain",
  which gives untransformed intervals, or "log", which gives
  log-transformed intervals.

- ci_level:

  Coverage probability of confidence intervals.

- transform:

  If set to `"log"`, `"logit"`, or `"odds"`, the standardized mean
  \\\theta(x)\\ is transformed into \\\psi(x)=\log\\\theta(x)\\\\,
  \\\psi(x)=\log\[\theta(x)/\\1-\theta(x)\\\]\\, or
  \\\psi(x)=\theta(x)/\\1-\theta(x)\\\\, respectively. If left
  unspecified, \\\psi(x)=\theta(x)\\.

- contrast:

  If set to `"difference"` or `"ratio"`, then \\\psi(x)-\psi(x_0)\\ or
  \\\psi(x) / \psi(x_0)\\ are constructed, where \\x_0\\ is a reference
  level specified by the `reference` argument. If not `NULL`, a doubly
  robust estimator of the standardized estimator is used.

- reference:

  If `contrast` is specified, the desired reference level.

- summary_fun:

  For internal use only. Do not change.

- ...:

  Unused.

## Value

None. Creates a plot as a side effect

## Examples

``` r
# see standardize_glm
```
