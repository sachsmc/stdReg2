# Implementing custom and new methods for standardization

``` r
library(stdReg2)
library(nnet)
```

## General implementation

The `stdReg2` package provides functionality for implementing regression
standardization for an arbitrary model. As long as the user can provide
a function to estimate the model, and another function to produce
predictions of the desired summary statistic for a given set of
confounders from the model, regression standardization can be applied to
produce estimates of a causal effect.

Here we demonstrate this by implementing a polytomous logistic
regression model for the outcome of activity in the nhefs data, which
takes three levels: (0) very active, (1) moderately active, and (2)
inactive. Inference can be done using the nonparametric bootstrap with
percentile confidence intervals (Tibshirani and Efron 1993) by setting
the argument `B` for the number of replicates, which we omit here for
brevity. Note that this method only applies for point exposures; other
methods are needed to produce valid effect estimates for time-varying
exposures, see Robins (1986) or the `gfoRmula` package (McGrath et al.
2020). Also note that the nonparametric bootstrap may not provide valid
inference if the outcome model is data-adaptive, e.g., based on machine
learning algorithms. In such situations alternative inference methods
may be required (Naimi, Mishler, and Kennedy 2021; Zivich and Breskin
2021).

``` r
nhefs_dat <- causaldata::nhefs_complete
# the target outcome model
# mfit <- multinom(active ~ qsmk + sex + race + age + I(age^2) + 
#               as.factor(education) + smokeintensity, data = nhefs_dat)

## here we predict the probability of being inactive (level 3)
predict_multinom <- function(...) {
  predict(..., type = "probs")[, 3]
}

std_custom <- standardize(fitter = "multinom", 
                          arguments = list(formula = active ~ qsmk + sex + 
                                             race + age + I(age^2) + 
               as.factor(education) + smokeintensity, trace = FALSE), 
               predict_fun = predict_multinom, 
               data = nhefs_dat, 
               values = list(qsmk = 0:1))
std_custom
#> Exposure:  qsmk 
#> Tables: 
#>  
#>   qsmk Estimate
#> 1    0     0.09
#> 2    1     0.11
```

To instruct `standardize` how to fit the target outcome model, we tell
it the name of the function to use in the `fitter` parameter
(`"multinom"` in this case) and provide the arguments to the fitter
function as a named list of `arguments`. The other required arguments
are the data, the values at which the causal effects are desired, and
the `predict_function`. The predict function takes as input the object
returned by fitter, and outputs a vector of predicted summary statistics
for each row of data. In this case we predict the probability of being
inactive, which is the third column of the result of predict on a
multinomial regression fit.

The output of the result is given above. We find the estimated
probability of being inactive is 9% in those who did not quit smoking,
while it is 11% in those who did quit smoking.

It is also possible to fit different models for each level of the
exposure, using the `standardize_level` function, as in the following
example:

``` r

std_custom2 <- standardize_level(fitter_list = list("multinom", "glm"),
                          arguments = list(list(formula = active ~ qsmk + sex + 
                                             race + age + I(age^2) + 
               as.factor(education) + smokeintensity, trace = FALSE), 
               list(formula = I(active == 2) ~ qsmk + sex + 
                                             race + age  + 
               as.factor(education) + smokeintensity, family = binomial)),
               predict_fun_list = list(predict_multinom, 
                                       \(...) predict.glm(..., type = "response")),
               data = nhefs_dat, 
               values = list(qsmk = 0:1))
std_custom2
#> Exposure:  qsmk 
#> Tables: 
#>  
#>   qsmk Estimate
#> 1    0   0.0894
#> 2    1   0.1117
```

Here we provide a list of fitters, arguments, and predict functions, one
for each of the desired levels of the exposure.

McGrath, Sean, Victoria Lin, Zilu Zhang, Lucia C. Petito, Roger W.
Logan, Miguel A. Hernán, and Jessica G. Young. 2020. “gfoRmula: An r
Package for Estimating the Effects of Sustained Treatment Strategies via
the Parametric g-Formula.” *Patterns* 1 (3): 100008.
https://doi.org/<https://doi.org/10.1016/j.patter.2020.100008>.

Naimi, Ashley I, Alan E Mishler, and Edward H Kennedy. 2021. “Challenges
in Obtaining Valid Causal Effect Estimates with Machine Learning
Algorithms.” *American Journal of Epidemiology* 192 (9): 1536–44.
<https://doi.org/10.1093/aje/kwab201>.

Robins, James. 1986. “A New Approach to Causal Inference in Mortality
Studies with a Sustained Exposure Period—Application to Control of the
Healthy Worker Survivor Effect.” *Mathematical Modelling* 7 (9-12):
1393–1512.

Tibshirani, Robert J, and Bradley Efron. 1993. “An Introduction to the
Bootstrap.” *Monographs on Statistics and Applied Probability* 57 (1):
1–436.

Zivich, Paul N., and Alexander Breskin. 2021. “Machine Learning for
Causal Inference: On the Use of Cross-Fit Estimators.” 32: 393–401.
<https://doi.org/10.1097/EDE.0000000000001332>.
