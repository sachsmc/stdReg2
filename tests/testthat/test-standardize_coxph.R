test_that("check estimates and standard errors against older package (coxph)", {
  require(survival)
  set.seed(8)
  n <- 300
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  time <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
  C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
  U <- pmin(time, C) # time at risk
  D <- as.numeric(time < C) # event indicator
  dd <- data.frame(Z, X, U, D)
  fit.std <- standardize_coxph(
    formula = Surv(U, D) ~ X + Z + X * Z,
    data = dd, values = list(X = seq(-1, 1, 0.5)), times = 3
  )

  expect_equal(unname(fit.std$res_contrast[[1]]$est_table[[1]]$Estimate), c(
    0.290086650197531, 0.119930794890419, 0.0812019338211797, 0.0796974035238777,
    0.0840306008325879
  ), tolerance = 1e-5)
  expect_equal(unname(fit.std$res_contrast[[1]]$est_table[[1]]$Std.Error), c(
    0.0550200346274478, 0.0365154454326768, 0.0265626695896798,
    0.0241162661987073, 0.0234315572499484
  ), tolerance = 1e-5)
})

