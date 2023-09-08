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

test_that("check estimates and standard errors standardize_glm (1)", {
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

