# Get regression standardized doubly-robust estimates from a glm

Get regression standardized doubly-robust estimates from a glm

## Usage

``` r
standardize_glm_dr(
  formula_outcome,
  formula_exposure,
  data,
  values,
  ci_level = 0.95,
  ci_type = "plain",
  contrasts = NULL,
  family_outcome = "gaussian",
  family_exposure = "binomial",
  reference = NULL,
  transforms = NULL
)
```

## Arguments

- formula_outcome:

  The formula which is used to fit the glm model for the outcome.

- formula_exposure:

  The formula which is used to fit the glm model for the exposure. If
  not `NULL`, a doubly robust estimator of the standardized estimator is
  used.

- data:

  The data.

- values:

  A named list or data.frame specifying the variables and values at
  which marginal means of the outcome will be estimated.

- ci_level:

  Coverage probability of confidence intervals.

- ci_type:

  A string, indicating the type of confidence intervals. Either "plain",
  which gives untransformed intervals, or "log", which gives
  log-transformed intervals.

- contrasts:

  A vector of contrasts in the following format: If set to
  `"difference"` or `"ratio"`, then \\\psi(x)-\psi(x_0)\\ or \\\psi(x) /
  \psi(x_0)\\ are constructed, where \\x_0\\ is a reference level
  specified by the `reference` argument. Has to be `NULL` if no
  references are specified.

- family_outcome:

  The family argument which is used to fit the glm model for the
  outcome.

- family_exposure:

  The family argument which is used to fit the glm model for the
  exposure.

- reference:

  A vector of reference levels in the following format: If `contrasts`
  is not `NULL`, the desired reference level(s). This must be a vector
  or list the same length as `contrasts`, and if not named, it is
  assumed that the order is as specified in contrasts.

- transforms:

  A vector of transforms in the following format: If set to `"log"`,
  `"logit"`, or `"odds"`, the standardized mean \\\theta(x)\\ is
  transformed into \\\psi(x)=\log\\\theta(x)\\\\,
  \\\psi(x)=\log\[\theta(x)/\\1-\theta(x)\\\]\\, or
  \\\psi(x)=\theta(x)/\\1-\theta(x)\\\\, respectively. If the vector is
  `NULL`, then \\\psi(x)=\theta(x)\\.

## Value

An object of class `std_glm`. Obtain numeric results in a data frame
with the
[tidy.std_glm](https://sachsmc.github.io/stdReg2/reference/tidy.std_glm.md)
function. This is a list with the following components:

- res_contrast:

  An unnamed list with one element for each of the requested contrasts.
  Each element is itself a list with the elements:

  estimates

  :   Estimated counterfactual means and standard errors for each
      exposure level

  covariance

  :   Estimated covariance matrix of counterfactual means

  fit_outcome

  :   The estimated regression model for the outcome

  fit_exposure

  :   The estimated exposure model

  exposure_names

  :   A character vector of the exposure variable names

  est_table

  :   Data.frame of the estimates of the contrast with inference

  transform

  :   The transform argument used for this contrast

  contrast

  :   The requested contrast type

  reference

  :   The reference level of the exposure

  ci_type

  :   Confidence interval type

  ci_level

  :   Confidence interval level

- res:

  A named list with the elements:

  estimates

  :   Estimated counterfactual means and standard errors for each
      exposure level

  covariance

  :   Estimated covariance matrix of counterfactual means

  fit_outcome

  :   The estimated regression model for the outcome

  fit_exposure

  :   The estimated exposure model

  exposure_names

  :   A character vector of the exposure variable names

## Details

`standardize_glm_dr` performs regression standardization in generalized
linear models, see e.g., documentation for `standardize_glm_dr`.
Specifically, this version uses a doubly robust estimator for
standardization, meaning inference is valid when either the outcome
regression or the exposure model is correctly specified and there is no
unmeasured confounding.

## References

Gabriel E.E., Sachs, M.C., Martinussen T., Waernbaum I., Goetghebeur E.,
Vansteelandt S., Sjölander A. (2024), Inverse probability of treatment
weighting with generalized linear outcome models for doubly robust
estimation. *Statistics in Medicine*, **43**(3):534–547.

## Examples

``` r
# doubly robust estimator
# needs to correctly specify either the outcome model or the exposure model
# for confounding
# NOTE: only works with binary exposures
x <- standardize_glm_dr(
  formula_outcome = bwt ~ smoker * (race + age + lwt) + I(age^2) + I(lwt^2),
  formula_exposure = smoker ~ race * age * lwt + I(age^2) + I(lwt^2),
  family_outcome = "gaussian",
  family_exposure = "binomial",
  data = clslowbwt,
  values = list(smoker = c(0, 1)), contrasts = "difference", reference = 0
)

set.seed(6)
n <- 100
Z <- rnorm(n)
X <- rbinom(n, 1, prob = (1 + exp(Z))^(-1))
Y <- rbinom(n, 1, prob = (1 + exp(as.numeric(X) + Z))^(-1))
dd <- data.frame(Z, X, Y)
x <- standardize_glm_dr(
  formula_outcome = Y ~ X * Z, formula_exposure = X ~ Z,
  family_outcome = "binomial",
  data = dd,
  values = list(X = 0:1), reference = 0,
  contrasts = c("difference"), transforms = c("odds")
)
```
