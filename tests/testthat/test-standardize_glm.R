test_that("bad names gives error", {
  set.seed(3)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean=Z)
  Y <- rnorm(n, mean=X+Z+0.1*X^2)
  dd <- data.frame(Z, X, Y)
  expect_error(standardize_glm(formula.outcome = Y~X+Z+I(X^2), fitbas, data=dd,values = list(X = 1, blue = 2)))
})

test_that("non-binary exposure gives error when using the doubly robust estimator", {
  set.seed(4)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean=Z)
  Y <- rnorm(n, mean=X+Z+0.1*X^2)
  dd <- data.frame(Z, X, Y)
  expect_error(standardize_glm(formula.outcome = Y~X+Z+I(X^2),formula.exposure = X~Z, fitbas, data=dd,values = list(X = c(-1,0,1))))
})

test_that("multiple exposure gives error when using the doubly robust estimator", {
  set.seed(5)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean=Z)
  X2 <- rnorm(n, mean=0.5*Z)
  Y <- rnorm(n, mean=X+Z+0.1*X^2)
  dd <- data.frame(Z, X,X2, Y)
  expect_error(standardize_glm(formula.outcome = Y~X+X2+Z+I(X^2),formula.exposure = X~Z, fitbas, data=dd,values = list(X = c(0,1), X2=c(0,1))))
})

test_that("check that it fails with p.population unspecified", {
  library(AF)
  data("singapore")
  Mi <- singapore$Age
  m <- mean(Mi)
  s <- sd(Mi)
  d <- 5
  expect_error(standardize_glm(formula.outcome = Oesophagealcancer~(Everhotbev+Age+ Dial+Samsu+Cigs)^2,
                  family.outcome=binomial, data=singapore,
                  values = list(Everhotbev = 0:1), clusterid = "Set",case.control = TRUE,
                  matched.density.cases = function(x) dnorm(x,m,s),
                  matched.density.controls = function(x) dnorm(x, m-d,s), matching.variable = Mi))
})

test_that("check estimates and standard errors standardize_glm (simple estimator)", {
  set.seed(6)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean=Z)
  Y <- rbinom(n, 1, prob=(1+exp(X+Z))^(-1))
  dd <- data.frame(Z, X, Y)
  # fit.stdReg1 <- glm(Y~X*Z, family=binomial,data=dd)
  # y <- stdGlm(fit.stdReg1,data=dd, X="X",x=0:1)
  # est.y <- summary(y)$est.table
  x <- standardize_glm(formula.outcome=Y~X*Z,family.outcome="binomial", data=dd, values = list(X = 0:1))
  expect_equal(x$estimates$estimates, c(0.519063874450474, 0.390531102199254), tolerance = 1e-5)
  expect_equal(x$estimates$se, c(0.0614996028608024, 0.0881636202817664), tolerance = 1e-5)
})

test_that("check estimates and standard errors standardize_glm (case-control estimator)",{
  library(AF)
  data("singapore")
  Mi <- singapore$Age
  m <- mean(Mi)
  s <- sd(Mi)
  d <- 5
  x<-suppressWarnings(standardize_glm(formula.outcome = Oesophagealcancer~(Everhotbev+Age+ Dial+Samsu+Cigs)^2,
                                      family.outcome=binomial, data=singapore,
                                      values = list(Everhotbev = 0:1), clusterid = "Set",case.control = TRUE,
                                      matched.density.cases = function(x) dnorm(x,m,s),
                                      matched.density.controls = function(x) dnorm(x, m-d,s), matching.variable = Mi,
                                      p.population = 19.3/100000))
  expect_equal(x$estimates$estimates, c(0.000127787218723152, 0.00057018810096904), tolerance = 1e-5)
  expect_equal(x$estimates$se, c(1.68549810412155e-05, 0.000218937483527179), tolerance = 1e-5)
})

test_that("check estimates and standard errors standardize_glm (case-control estimator)",{
  library(AF)
  data <- AF::clslowbwt
  x<-standardize_glm(formula.outcome = bwt ~ smoker * (race + age + lwt) + I(age^2) + I(lwt^2),
                     formula.exposure = smoker ~ race * age * lwt + I(age^2) + I(lwt^2),
                     family.outcome = gaussian,
                     family.exposure = binomial,
                     data = data,
                     values = list(smoker = c(0,1)))
  z<-summary(x,contrast="difference",reference = 0)

  expect_equal(as.numeric(z$est.table$Estimate)[2], -223.6736, tolerance = 1e-5)
  expect_equal(as.numeric(z$est.table$`lower 0.95`)[2], -446.0778, tolerance = 1e-5)
  expect_equal(as.numeric(z$est.table$`upper 0.95`)[2], -1.269345, tolerance = 1e-5)
})

