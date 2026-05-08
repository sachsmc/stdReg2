# Get standardized estimates using the g-formula with a custom model

Get standardized estimates using the g-formula with a custom model

## Usage

``` r
standardize(
  fitter,
  arguments,
  predict_fun,
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

- fitter:

  The function to call to fit the data.

- arguments:

  The arguments to be used in the fitter function as a `list`.

- predict_fun:

  The function used to predict the means/probabilities for a new data
  set on the response level. For survival data, this should be a matrix
  where each column is the time, and each row the data.

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

Let \\Y\\, \\X\\, and \\Z\\ be the outcome, the exposure, and a vector
of covariates, respectively. `standardize` uses a model to estimate the
standardized mean \\\theta(x)=E\\E(Y\|X=x,Z)\\\\, where \\x\\ is a
specific value of \\X\\, and the outer expectation is over the marginal
distribution of \\Z\\. With survival data, \\Y=I(T \> t)\\, and a vector
of different time points `times` (\\t\\) can be given, where \\T\\ is
the uncensored survival time.

Note that the nonparametric bootstrap may not provide valid inference if
the outcome model is data-adaptive, e.g., based on machine learning
algorithms. In such situations alternative inference methods may be
required.

## References

Rothman K.J., Greenland S., Lash T.L. (2008). *Modern Epidemiology*, 3rd
edition. Lippincott, Williams & Wilkins.

Sjölander A. (2016). Regression standardization with the R-package
stdReg. *European Journal of Epidemiology* **31**(6), 563-574.

Sjölander A. (2016). Estimation of causal effect measures with the
R-package stdReg. *European Journal of Epidemiology* **33**(9), 847-858.

## Examples

``` r
set.seed(6)
n <- 100
Z <- rnorm(n)
X <- rnorm(n, mean = Z)
Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
dd <- data.frame(Z, X, Y)
prob_predict.glm <- function(...) predict.glm(..., type = "response")

x <- standardize(
  fitter = "glm",
  arguments = list(
    formula = Y ~ X * Z,
    family = "binomial"
  ),
  predict_fun = prob_predict.glm,
  data = dd,
  values = list(X = seq(-1, 1, 0.1)),
  B = 100,
  reference = 0,
  contrasts = "difference"
)
#>   |                                                          |                                                  |   0%  |                                                          |=                                                 |   1%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   3%  |                                                          |==                                                |   4%  |                                                          |===                                               |   5%  |                                                          |===                                               |   6%  |                                                          |====                                              |   7%  |                                                          |====                                              |   8%  |                                                          |=====                                             |   9%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  11%  |                                                          |======                                            |  12%  |                                                          |=======                                           |  13%  |                                                          |=======                                           |  14%  |                                                          |========                                          |  15%  |                                                          |========                                          |  16%  |                                                          |=========                                         |  17%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  19%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  21%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  23%  |                                                          |============                                      |  24%  |                                                          |=============                                     |  25%  |                                                          |=============                                     |  26%  |                                                          |==============                                    |  27%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  29%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  31%  |                                                          |================                                  |  32%  |                                                          |=================                                 |  33%  |                                                          |=================                                 |  34%  |                                                          |==================                                |  35%  |                                                          |==================                                |  36%  |                                                          |===================                               |  37%  |                                                          |===================                               |  38%  |                                                          |====================                              |  39%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  41%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  43%  |                                                          |======================                            |  44%  |                                                          |=======================                           |  45%  |                                                          |=======================                           |  46%  |                                                          |========================                          |  47%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  49%  |                                                          |=========================                         |  51%  |                                                          |==========================                        |  52%  |                                                          |==========================                        |  53%  |                                                          |===========================                       |  54%  |                                                          |===========================                       |  55%  |                                                          |============================                      |  56%  |                                                          |============================                      |  57%  |                                                          |=============================                     |  58%  |                                                          |=============================                     |  59%  |                                                          |==============================                    |  60%  |                                                          |==============================                    |  61%  |                                                          |===============================                   |  62%  |                                                          |===============================                   |  63%  |                                                          |================================                  |  64%  |                                                          |================================                  |  65%  |                                                          |=================================                 |  66%  |                                                          |=================================                 |  67%  |                                                          |==================================                |  68%  |                                                          |==================================                |  69%  |                                                          |===================================               |  70%  |                                                          |===================================               |  71%  |                                                          |====================================              |  72%  |                                                          |====================================              |  73%  |                                                          |=====================================             |  74%  |                                                          |=====================================             |  75%  |                                                          |======================================            |  76%  |                                                          |======================================            |  77%  |                                                          |=======================================           |  78%  |                                                          |=======================================           |  79%  |                                                          |========================================          |  80%  |                                                          |========================================          |  81%  |                                                          |=========================================         |  82%  |                                                          |=========================================         |  83%  |                                                          |==========================================        |  84%  |                                                          |==========================================        |  85%  |                                                          |===========================================       |  86%  |                                                          |===========================================       |  87%  |                                                          |============================================      |  88%  |                                                          |============================================      |  89%  |                                                          |=============================================     |  90%  |                                                          |=============================================     |  91%  |                                                          |==============================================    |  92%  |                                                          |==============================================    |  93%  |                                                          |===============================================   |  94%  |                                                          |===============================================   |  95%  |                                                          |================================================  |  96%  |                                                          |================================================  |  97%  |                                                          |================================================= |  98%  |                                                          |================================================= |  99%  |                                                          |==================================================| 100%
x
#> Number of bootstraps:  100 
#> Confidence intervals are based on percentile bootstrap confidence intervals 
#> 
#> Exposure:  X 
#> Tables: 
#>  
#>       X Estimate lower.0.95 upper.0.95
#> 1  -1.0    0.689      0.537      0.882
#> 2  -0.9    0.671      0.526      0.865
#> 3  -0.8    0.653      0.521      0.845
#> 4  -0.7    0.635      0.515      0.822
#> 5  -0.6    0.617      0.508      0.797
#> 6  -0.5    0.600      0.493      0.770
#> 7  -0.4    0.583      0.476      0.745
#> 8  -0.3    0.566      0.466      0.720
#> 9  -0.2    0.550      0.450      0.699
#> 10 -0.1    0.534      0.429      0.679
#> 11  0.0    0.519      0.409      0.658
#> 12  0.1    0.504      0.389      0.638
#> 13  0.2    0.490      0.369      0.621
#> 14  0.3    0.476      0.348      0.608
#> 15  0.4    0.462      0.322      0.595
#> 16  0.5    0.449      0.298      0.587
#> 17  0.6    0.437      0.275      0.583
#> 18  0.7    0.425      0.254      0.577
#> 19  0.8    0.413      0.235      0.571
#> 20  0.9    0.401      0.222      0.561
#> 21  1.0    0.391      0.207      0.551
#> 
#> Reference level:  = 0 
#> Contrast:  difference 
#>       X Estimate lower.0.95 upper.0.95
#> 1  -1.0   0.1697    0.04037    0.30879
#> 2  -0.9   0.1516    0.03562    0.28505
#> 3  -0.8   0.1335    0.03104    0.25348
#> 4  -0.7   0.1156    0.02663    0.22134
#> 5  -0.6   0.0980    0.02238    0.19201
#> 6  -0.5   0.0807    0.01829    0.16211
#> 7  -0.4   0.0637    0.01435    0.13054
#> 8  -0.3   0.0472    0.01056    0.09810
#> 9  -0.2   0.0310    0.00691    0.06525
#> 10 -0.1   0.0153    0.00339    0.03243
#> 11  0.0   0.0000    0.00000    0.00000
#> 12  0.1  -0.0148   -0.03177   -0.00326
#> 13  0.2  -0.0292   -0.06286   -0.00641
#> 14  0.3  -0.0432   -0.09365   -0.00944
#> 15  0.4  -0.0566   -0.12370   -0.01237
#> 16  0.5  -0.0697   -0.15192   -0.01519
#> 17  0.6  -0.0823   -0.17663   -0.01791
#> 18  0.7  -0.0945   -0.19954   -0.02054
#> 19  0.8  -0.1062   -0.22074   -0.02308
#> 20  0.9  -0.1176   -0.24033   -0.02553
#> 21  1.0  -0.1285   -0.25891   -0.02790
#> 

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
X <- rnorm(n, mean = Z)
T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
U <- pmin(T, C) # time at risk
D <- as.numeric(T < C) # event indicator
dd <- data.frame(Z, X, U, D)
x <- standardize(
fitter = "coxph",
  arguments = list(
    formula = Surv(U, D) ~ X + Z + X * Z,
    method = "breslow",
    x = TRUE,
    y = TRUE
  ),
  predict_fun = prob_predict.coxph,
  data = dd,
  times = 1:5,
  values = list(X = c(-1, 0, 1)),
  B = 100,
  reference = 0,
  contrasts = "difference"
)
#>   |                                                          |                                                  |   0%  |                                                          |=                                                 |   1%  |                                                          |=                                                 |   2%  |                                                          |==                                                |   3%  |                                                          |==                                                |   4%  |                                                          |===                                               |   5%  |                                                          |===                                               |   6%  |                                                          |====                                              |   7%  |                                                          |====                                              |   8%  |                                                          |=====                                             |   9%  |                                                          |=====                                             |  10%  |                                                          |======                                            |  11%  |                                                          |======                                            |  12%  |                                                          |=======                                           |  13%  |                                                          |=======                                           |  14%  |                                                          |========                                          |  15%  |                                                          |========                                          |  16%  |                                                          |=========                                         |  17%  |                                                          |=========                                         |  18%  |                                                          |==========                                        |  19%  |                                                          |==========                                        |  20%  |                                                          |===========                                       |  21%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  23%  |                                                          |============                                      |  24%  |                                                          |=============                                     |  25%  |                                                          |=============                                     |  26%  |                                                          |==============                                    |  27%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  29%  |                                                          |===============                                   |  30%  |                                                          |================                                  |  31%  |                                                          |================                                  |  32%  |                                                          |=================                                 |  33%  |                                                          |=================                                 |  34%  |                                                          |==================                                |  35%  |                                                          |==================                                |  36%  |                                                          |===================                               |  37%  |                                                          |===================                               |  38%  |                                                          |====================                              |  39%  |                                                          |====================                              |  40%  |                                                          |=====================                             |  41%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  43%  |                                                          |======================                            |  44%  |                                                          |=======================                           |  45%  |                                                          |=======================                           |  46%  |                                                          |========================                          |  47%  |                                                          |========================                          |  48%  |                                                          |=========================                         |  49%  |                                                          |=========================                         |  51%  |                                                          |==========================                        |  52%  |                                                          |==========================                        |  53%  |                                                          |===========================                       |  54%  |                                                          |===========================                       |  55%  |                                                          |============================                      |  56%  |                                                          |============================                      |  57%  |                                                          |=============================                     |  58%  |                                                          |=============================                     |  59%  |                                                          |==============================                    |  60%  |                                                          |==============================                    |  61%  |                                                          |===============================                   |  62%  |                                                          |===============================                   |  63%  |                                                          |================================                  |  64%  |                                                          |================================                  |  65%  |                                                          |=================================                 |  66%  |                                                          |=================================                 |  67%  |                                                          |==================================                |  68%  |                                                          |==================================                |  69%  |                                                          |===================================               |  70%  |                                                          |===================================               |  71%  |                                                          |====================================              |  72%  |                                                          |====================================              |  73%  |                                                          |=====================================             |  74%  |                                                          |=====================================             |  75%  |                                                          |======================================            |  76%  |                                                          |======================================            |  77%  |                                                          |=======================================           |  78%  |                                                          |=======================================           |  79%  |                                                          |========================================          |  80%  |                                                          |========================================          |  81%  |                                                          |=========================================         |  82%  |                                                          |=========================================         |  83%  |                                                          |==========================================        |  84%  |                                                          |==========================================        |  85%  |                                                          |===========================================       |  86%  |                                                          |===========================================       |  87%  |                                                          |============================================      |  88%  |                                                          |============================================      |  89%  |                                                          |=============================================     |  90%  |                                                          |=============================================     |  91%  |                                                          |==============================================    |  92%  |                                                          |==============================================    |  93%  |                                                          |===============================================   |  94%  |                                                          |===============================================   |  95%  |                                                          |================================================  |  96%  |                                                          |================================================  |  97%  |                                                          |================================================= |  98%  |                                                          |================================================= |  99%  |                                                          |==================================================| 100%
x
#> Number of bootstraps:  100 
#> Confidence intervals are based on percentile bootstrap confidence intervals 
#> 
#> Exposure:  X 
#> Tables: 
#>  
#> Time:  1 
#>    X  estimate     lower     upper
#> 1 -1 0.6691328 0.6078557 0.7171407
#> 2  0 0.3632410 0.3011138 0.4060458
#> 3  1 0.2466151 0.1878811 0.2902875
#> 
#> Time:  2 
#>    X  estimate     lower     upper
#> 1 -1 0.4676722 0.3847559 0.5519054
#> 2  0 0.2135942 0.1659279 0.2611208
#> 3  1 0.1656470 0.1168408 0.2079257
#> 
#> Time:  3 
#>    X  estimate      lower     upper
#> 1 -1 0.3505275 0.26905828 0.4459009
#> 2  0 0.1520739 0.11902194 0.2044295
#> 3  1 0.1315275 0.09380255 0.1761028
#> 
#> Time:  4 
#>    X  estimate      lower     upper
#> 1 -1 0.2546911 0.12386946 0.3508293
#> 2  0 0.1101689 0.06701637 0.1505347
#> 3  1 0.1069771 0.06241746 0.1430167
#> 
#> Time:  5 
#>    X   estimate      lower     upper
#> 1 -1 0.16309819 0.07709911 0.2613375
#> 2  0 0.07480425 0.03816440 0.1127717
#> 3  1 0.08456621 0.04640860 0.1207358
#> 
#> 
#> Reference level:  = 0 
#> Contrast:  difference 
#> Time:  1 
#>    X   estimate      lower       upper
#> 1 -1  0.3058918  0.2515877  0.35267464
#> 2  0  0.0000000  0.0000000  0.00000000
#> 3  1 -0.1166259 -0.1534746 -0.08540127
#> 
#> Time:  2 
#>    X    estimate       lower       upper
#> 1 -1  0.25407801  0.18502799  0.33076959
#> 2  0  0.00000000  0.00000000  0.00000000
#> 3  1 -0.04794719 -0.07480939 -0.02350433
#> 
#> Time:  3 
#>    X    estimate       lower       upper
#> 1 -1  0.19845364  0.12283854 0.284032243
#> 2  0  0.00000000  0.00000000 0.000000000
#> 3  1 -0.02054638 -0.04719901 0.001297897
#> 
#> Time:  4 
#>    X     estimate       lower      upper
#> 1 -1  0.144522281  0.04976681 0.22983546
#> 2  0  0.000000000  0.00000000 0.00000000
#> 3  1 -0.003191768 -0.02345454 0.01911037
#> 
#> Time:  5 
#>    X    estimate        lower      upper
#> 1 -1 0.088293936  0.018981386 0.17857030
#> 2  0 0.000000000  0.000000000 0.00000000
#> 3  1 0.009761965 -0.009490497 0.02585953
#> 
#> 
```
