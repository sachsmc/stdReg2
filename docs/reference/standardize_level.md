# Get standardized estimates using the g-formula with and separate models for each exposure level in the data

Get standardized estimates using the g-formula with and separate models
for each exposure level in the data

## Usage

``` r
standardize_level(
  fitter_list,
  arguments,
  predict_fun_list,
  data,
  values,
  B = NULL,
  ci_level = 0.95,
  contrasts = NULL,
  reference = NULL,
  seed = NULL,
  times = NULL,
  transforms = NULL,
  progressbar = TRUE
)
```

## Arguments

- fitter_list:

  The function to call to fit the data (as a list).

- arguments:

  The arguments to be used in the fitter function as a `list`.

- predict_fun_list:

  The function used to predict the means/probabilities for a new data
  set on the response level. For survival data, this should be a matrix
  where each column is the time, and each row the data (as a list).

- data:

  The data.

- values:

  A named list or data.frame specifying the variables and values at
  which marginal means of the outcome will be estimated.

- B:

  Number of nonparametric bootstrap resamples. Default is `NULL` (no
  bootstrap).

- ci_level:

  Coverage probability of confidence intervals.

- contrasts:

  A vector of contrasts in the following format: If set to
  `"difference"` or `"ratio"`, then \\\psi(x)-\psi(x_0)\\ or \\\psi(x) /
  \psi(x_0)\\ are constructed, where \\x_0\\ is a reference level
  specified by the `reference` argument. Has to be `NULL` if no
  references are specified.

- reference:

  A vector of reference levels in the following format: If `contrasts`
  is not `NULL`, the desired reference level(s). This must be a vector
  or list the same length as `contrasts`, and if not named, it is
  assumed that the order is as specified in contrasts.

- seed:

  The seed to use with the nonparametric bootstrap.

- times:

  For use with survival data. Set to `NULL` otherwise.

- transforms:

  A vector of transforms in the following format: If set to `"log"`,
  `"logit"`, or `"odds"`, the standardized mean \\\theta(x)\\ is
  transformed into \\\psi(x)=\log\\\theta(x)\\\\,
  \\\psi(x)=\log\[\theta(x)/\\1-\theta(x)\\\]\\, or
  \\\psi(x)=\theta(x)/\\1-\theta(x)\\\\, respectively. If the vector is
  `NULL`, then \\\psi(x)=\theta(x)\\.

- progressbar:

  Logical, if TRUE will print bootstrapping progress to the console

## Value

An object of class `std_custom`. Obtain numeric results using
[tidy.std_custom](https://sachsmc.github.io/stdReg2/reference/tidy.std_custom.md).
This is a list with the following components:

- res_contrast:

  An unnamed list with one element for each of the requested contrasts.
  Each element is itself a list with the elements:

  B

  :   The number of bootstrap replicates

  estimates

  :   Estimated counterfactual means and standard errors for each
      exposure level

  fit_outcome

  :   The estimated regression model for the outcome

  estimates_boot

  :   A list of estimates, one for each bootstrap resample

  exposure_names

  :   A character vector of the exposure variable names

  times

  :   The vector of times at which the calculation is done, if relevant

  est_table

  :   Data.frame of the estimates of the contrast with inference

  transform

  :   The transform argument used for this contrast

  contrast

  :   The requested contrast type

  reference

  :   The reference level of the exposure

  ci_level

  :   Confidence interval level

- res:

  A named list with the elements:

  B

  :   The number of bootstrap replicates

  estimates

  :   Estimated counterfactual means and standard errors for each
      exposure level

  fit_outcome

  :   The estimated regression model for the outcome

  estimates_boot

  :   A list of estimates, one for each bootstrap resample

  exposure_names

  :   A character vector of the exposure variable names

  times

  :   The vector of times at which the calculation is done, if relevant

## Details

See `standardize`. The difference is here that different models can be
fitted for each value of `x` in `values`.

## References

Rothman K.J., Greenland S., Lash T.L. (2008). *Modern Epidemiology*, 3rd
edition. Lippincott, Williams & Wilkins.

Sjölander A. (2016). Regression standardization with the R-package
stdReg. *European Journal of Epidemiology* **31**(6), 563-574.

Sjölander A. (2016). Estimation of causal effect measures with the
R-package stdReg. *European Journal of Epidemiology* **33**(9), 847-858.

## Examples

``` r

require(survival)
prob_predict.coxph <- function(object, newdata, times) {
  fit.detail <- suppressWarnings(basehaz(object))
  cum.haz <- fit.detail$hazard[sapply(times, function(x) max(which(fit.detail$time <= x)))]
  predX <- predict(object = object, newdata = newdata, type = "risk")
  res <- matrix(NA, ncol = length(times), nrow = length(predX))
  for (ti in seq_len(length(times))) {
    res[, ti] <- exp(-predX * cum.haz[ti])
  }
  res
}

set.seed(68)
n <- 500
Z <- rnorm(n)
X <- rbinom(n, 1, prob = 0.5)
T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
U <- pmin(T, C) # time at risk
D <- as.numeric(T < C) # event indicator
dd <- data.frame(Z, X, U, D)
x <- standardize_level(
  fitter_list = list("coxph", "coxph"),
  arguments = list(
    list(
      formula = Surv(U, D) ~ X + Z + X * Z,
      method = "breslow",
      x = TRUE,
      y = TRUE
    ),
    list(
      formula = Surv(U, D) ~ X,
      method = "breslow",
      x = TRUE,
      y = TRUE
    )
  ),
  predict_fun_list = list(prob_predict.coxph, prob_predict.coxph),
  data = dd,
  times = seq(1, 5, 0.1),
  values = list(X = c(0, 1)),
  B = 100,
  reference = 0,
  contrasts = "difference"
)
#>   |                                                          |                                                  |   0%  |                                                          |=                                                 |   1%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   3%  |                                                          |==                                                |   4%  |                                                          |===                                               |   5%  |                                                          |===                                               |   6%  |                                                          |====                                              |   7%  |                                                          |====                                              |   8%  |                                                          |=====                                             |   9%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  11%  |                                                          |======                                            |  12%  |                                                          |=======                                           |  13%  |                                                          |=======                                           |  14%  |                                                          |========                                          |  15%  |                                                          |========                                          |  16%  |                                                          |=========                                         |  17%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  19%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  21%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  23%  |                                                          |============                                      |  24%  |                                                          |=============                                     |  25%  |                                                          |=============                                     |  26%  |                                                          |==============                                    |  27%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  29%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  31%  |                                                          |================                                  |  32%  |                                                          |=================                                 |  33%  |                                                          |=================                                 |  34%  |                                                          |==================                                |  35%  |                                                          |==================                                |  36%  |                                                          |===================                               |  37%  |                                                          |===================                               |  38%  |                                                          |====================                              |  39%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  41%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  43%  |                                                          |======================                            |  44%  |                                                          |=======================                           |  45%  |                                                          |=======================                           |  46%  |                                                          |========================                          |  47%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  49%  |                                                          |=========================                         |  51%  |                                                          |==========================                        |  52%  |                                                          |==========================                        |  53%  |                                                          |===========================                       |  54%  |                                                          |===========================                       |  55%  |                                                          |============================                      |  56%  |                                                          |============================                      |  57%  |                                                          |=============================                     |  58%  |                                                          |=============================                     |  59%  |                                                          |==============================                    |  60%  |                                                          |==============================                    |  61%  |                                                          |===============================                   |  62%  |                                                          |===============================                   |  63%  |                                                          |================================                  |  64%  |                                                          |================================                  |  65%  |                                                          |=================================                 |  66%  |                                                          |=================================                 |  67%  |                                                          |==================================                |  68%  |                                                          |==================================                |  69%  |                                                          |===================================               |  70%  |                                                          |===================================               |  71%  |                                                          |====================================              |  72%  |                                                          |====================================              |  73%  |                                                          |=====================================             |  74%  |                                                          |=====================================             |  75%  |                                                          |======================================            |  76%  |                                                          |======================================            |  77%  |                                                          |=======================================           |  78%  |                                                          |=======================================           |  79%  |                                                          |========================================          |  80%  |                                                          |========================================          |  81%  |                                                          |=========================================         |  82%  |                                                          |=========================================         |  83%  |                                                          |==========================================        |  84%  |                                                          |==========================================        |  85%  |                                                          |===========================================       |  86%  |                                                          |===========================================       |  87%  |                                                          |============================================      |  88%  |                                                          |============================================      |  89%  |                                                          |=============================================     |  90%  |                                                          |=============================================     |  91%  |                                                          |==============================================    |  92%  |                                                          |==============================================    |  93%  |                                                          |===============================================   |  94%  |                                                          |===============================================   |  95%  |                                                          |================================================  |  96%  |                                                          |================================================  |  97%  |                                                          |================================================= |  98%  |                                                          |================================================= |  99%  |                                                          |==================================================| 100%
print(x)
#> Number of bootstraps:  100 
#> Confidence intervals are based on percentile bootstrap confidence intervals 
#> 
#> Exposure:  X 
#> Tables: 
#>  
#> Time:  1 
#>   X  estimate     lower     upper
#> 1 0 0.3426359 0.2685875 0.4126004
#> 2 1 0.3306066 0.2719930 0.3968803
#> 
#> Time:  1.1 
#>   X  estimate     lower     upper
#> 1 0 0.3193937 0.2478319 0.3800825
#> 2 1 0.3140225 0.2554192 0.3785019
#> 
#> Time:  1.2 
#>   X  estimate     lower     upper
#> 1 0 0.2910216 0.2229411 0.3537408
#> 2 1 0.2929892 0.2427978 0.3563471
#> 
#> Time:  1.3 
#>   X  estimate     lower     upper
#> 1 0 0.2784148 0.2174233 0.3452300
#> 2 1 0.2844500 0.2347156 0.3465131
#> 
#> Time:  1.4 
#>   X  estimate     lower     upper
#> 1 0 0.2722029 0.2119961 0.3390574
#> 2 1 0.2801398 0.2333808 0.3414379
#> 
#> Time:  1.5 
#>   X  estimate     lower     upper
#> 1 0 0.2537637 0.1902978 0.3208143
#> 2 1 0.2670585 0.2136579 0.3222557
#> 
#> Time:  1.6 
#>   X  estimate     lower     upper
#> 1 0 0.2473441 0.1836796 0.3202737
#> 2 1 0.2625318 0.2112083 0.3187268
#> 
#> Time:  1.7 
#>   X  estimate     lower     upper
#> 1 0 0.2284202 0.1650471 0.2935825
#> 2 1 0.2487069 0.1926763 0.3090800
#> 
#> Time:  1.8 
#>   X  estimate     lower     upper
#> 1 0 0.2284202 0.1650471 0.2935825
#> 2 1 0.2487069 0.1926763 0.3090800
#> 
#> Time:  1.9 
#>   X  estimate     lower     upper
#> 1 0 0.2155678 0.1545570 0.2795913
#> 2 1 0.2385031 0.1879035 0.2981300
#> 
#> Time:  2 
#>   X  estimate     lower     upper
#> 1 0 0.2083187 0.1540708 0.2697798
#> 2 1 0.2328579 0.1835953 0.2894313
#> 
#> Time:  2.1 
#>   X  estimate     lower     upper
#> 1 0 0.2010004 0.1422739 0.2664994
#> 2 1 0.2271435 0.1820213 0.2827236
#> 
#> Time:  2.2 
#>   X  estimate     lower     upper
#> 1 0 0.2010004 0.1422739 0.2664994
#> 2 1 0.2271435 0.1820213 0.2827236
#> 
#> Time:  2.3 
#>   X  estimate     lower     upper
#> 1 0 0.1931929 0.1335466 0.2582027
#> 2 1 0.2210044 0.1767496 0.2766781
#> 
#> Time:  2.4 
#>   X  estimate     lower     upper
#> 1 0 0.1931929 0.1335466 0.2582027
#> 2 1 0.2210044 0.1767496 0.2766781
#> 
#> Time:  2.5 
#>   X  estimate     lower     upper
#> 1 0 0.1648695 0.1077254 0.2294345
#> 2 1 0.1993885 0.1465611 0.2613009
#> 
#> Time:  2.6 
#>   X  estimate     lower     upper
#> 1 0 0.1557581 0.1009572 0.2232631
#> 2 1 0.1920041 0.1396859 0.2496035
#> 
#> Time:  2.7 
#>   X  estimate     lower     upper
#> 1 0 0.1557581 0.1009572 0.2232631
#> 2 1 0.1920041 0.1396859 0.2496035
#> 
#> Time:  2.8 
#>   X  estimate     lower     upper
#> 1 0 0.1557581 0.1009572 0.2232631
#> 2 1 0.1920041 0.1396859 0.2496035
#> 
#> Time:  2.9 
#>   X  estimate     lower     upper
#> 1 0 0.1557581 0.1009572 0.2232631
#> 2 1 0.1920041 0.1396859 0.2496035
#> 
#> Time:  3 
#>   X  estimate      lower     upper
#> 1 0 0.1455376 0.09859477 0.2099133
#> 2 1 0.1838344 0.13550757 0.2463232
#> 
#> Time:  3.1 
#>   X  estimate      lower     upper
#> 1 0 0.1344841 0.08736479 0.1989949
#> 2 1 0.1752849 0.13259212 0.2324468
#> 
#> Time:  3.2 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  3.3 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  3.4 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  3.5 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  3.6 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  3.7 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  3.8 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  3.9 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.1 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.2 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.3 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.4 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.5 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.6 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.7 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.8 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  4.9 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> Time:  5 
#>   X  estimate      lower     upper
#> 1 0 0.1239428 0.08002961 0.1812128
#> 2 1 0.1662973 0.12094475 0.2177145
#> 
#> 
#> Reference level:  = 0 
#> Contrast:  difference 
#> Time:  1 
#>   X    estimate       lower    upper
#> 1 0  0.00000000  0.00000000 0.000000
#> 2 1 -0.01202936 -0.08167224 0.063079
#> 
#> Time:  1.1 
#>   X     estimate       lower      upper
#> 1 0  0.000000000  0.00000000 0.00000000
#> 2 1 -0.005371227 -0.07540159 0.07050428
#> 
#> Time:  1.2 
#>   X    estimate       lower      upper
#> 1 0 0.000000000  0.00000000 0.00000000
#> 2 1 0.001967654 -0.06990685 0.07668053
#> 
#> Time:  1.3 
#>   X    estimate      lower     upper
#> 1 0 0.000000000  0.0000000 0.0000000
#> 2 1 0.006035151 -0.0652183 0.0817649
#> 
#> Time:  1.4 
#>   X    estimate       lower      upper
#> 1 0 0.000000000  0.00000000 0.00000000
#> 2 1 0.007936927 -0.06470686 0.08250201
#> 
#> Time:  1.5 
#>   X   estimate       lower    upper
#> 1 0 0.00000000  0.00000000 0.000000
#> 2 1 0.01329481 -0.05696887 0.086919
#> 
#> Time:  1.6 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.01518772 -0.05696887 0.08771999
#> 
#> Time:  1.7 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.02028668 -0.04978553 0.09089055
#> 
#> Time:  1.8 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.02028668 -0.04978553 0.09089055
#> 
#> Time:  1.9 
#>   X   estimate       lower     upper
#> 1 0 0.00000000  0.00000000 0.0000000
#> 2 1 0.02293534 -0.04767694 0.0918496
#> 
#> Time:  2 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.02453919 -0.04613972 0.09298423
#> 
#> Time:  2.1 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.02614307 -0.04472088 0.09315902
#> 
#> Time:  2.2 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.02614307 -0.04472088 0.09315902
#> 
#> Time:  2.3 
#>   X  estimate       lower      upper
#> 1 0 0.0000000  0.00000000 0.00000000
#> 2 1 0.0278115 -0.04234635 0.09415061
#> 
#> Time:  2.4 
#>   X  estimate       lower      upper
#> 1 0 0.0000000  0.00000000 0.00000000
#> 2 1 0.0278115 -0.04234635 0.09415061
#> 
#> Time:  2.5 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.03451898 -0.03573643 0.09685651
#> 
#> Time:  2.6 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.03624594 -0.03573643 0.09685959
#> 
#> Time:  2.7 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.03624594 -0.03573643 0.09685959
#> 
#> Time:  2.8 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.03624594 -0.03573643 0.09685959
#> 
#> Time:  2.9 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.03624594 -0.03573643 0.09685959
#> 
#> Time:  3 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.03829678 -0.03307009 0.09632176
#> 
#> Time:  3.1 
#>   X  estimate       lower      upper
#> 1 0 0.0000000  0.00000000 0.00000000
#> 2 1 0.0408008 -0.02680682 0.09743773
#> 
#> Time:  3.2 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  3.3 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  3.4 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  3.5 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  3.6 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  3.7 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  3.8 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  3.9 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.1 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.2 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.3 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.4 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.5 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.6 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.7 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.8 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  4.9 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> Time:  5 
#>   X   estimate       lower      upper
#> 1 0 0.00000000  0.00000000 0.00000000
#> 2 1 0.04235446 -0.02413252 0.09619845
#> 
#> 
```
