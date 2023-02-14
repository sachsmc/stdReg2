test_that("bad names gives error", {
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean=Z)
  Y <- rnorm(n, mean=X+Z+0.1*X^2)
  dd <- data.frame(Z, X, Y)
  fitbas <- glm(formula=Y~X+Z+I(X^2), data=dd)

  expect_error(standardize.glm(fitbas, values = list(X = 1, blue = 2)))
})
