test_that("bad names gives error", {
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean=Z)
  Y <- rnorm(n, mean=X+Z+0.1*X^2)
  dd <- data.frame(Z, X, Y)
  expect_error(standardize_glm(formula.outcome = Y~X+Z+I(X^2), fitbas, data=dd,values = list(X = 1, blue = 2)))
})
