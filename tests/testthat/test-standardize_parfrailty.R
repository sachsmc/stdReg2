test_that("check estimates and standard errors against older package (parfrailty)", {
  require(survival)
  # simulate data
  set.seed(6)
  n <- 300
  m <- 3
  alpha <- 1.5
  eta <- 1
  phi <- 0.5
  beta <- 1
  id <- rep(1:n, each = m)
  U <- rep(rgamma(n, shape = 1 / phi, scale = phi), each = m)
  X <- rnorm(n * m)
  # reparametrize scale as in rweibull function
  weibull.scale <- alpha / (U * exp(beta * X))^(1 / eta)
  T <- rweibull(n * m, shape = eta, scale = weibull.scale)

  # right censoring
  C <- runif(n * m, 0, 10)
  D <- as.numeric(T < C)
  T <- pmin(T, C)

  # strong left-truncation
  L <- runif(n * m, 0, 2)
  incl <- T > L
  incl <- ave(x = incl, id, FUN = sum) == m
  dd <- data.frame(L, T, D, X, id)
  dd <- dd[incl, ]

  fit <- parfrailty(formula = Surv(L, T, D) ~ X, data = dd, clusterid = "id")
  fit.std <- stdParfrailty(fit = fit, data = dd, X = "X", x = seq(-1, 1, 0.5), t = 3, clusterid = "id")
  x <- summary(fit.std)

  expect_equal(unname(x$est.table[[1]][, 1]), c(
    0.568209706629598, 0.446427416102044, 0.322741943035805, 0.211677902336095,
    0.124866721293621
  ), tolerance = 1e-5)
  expect_equal(unname(x$est.table[[1]][, 2]), c(
    0.069294417658677, 0.0700174182939671, 0.0685777998185612,
    0.0632225540159503, 0.0525973717447175
  ), tolerance = 1e-5)
})