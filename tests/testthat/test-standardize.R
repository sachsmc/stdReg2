test_that("standardize generates output with non-survival data (no bootstrap)", {
  set.seed(6)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
  dd <- data.frame(Z, X, Y)
  prob_predict.glm <- function(...) predict.glm(..., type = "response")

  x <- standardize(
    arguments = list(
      formula = Y ~ X * Z,
      family = "binomial"
    ),
    predict_fun = prob_predict.glm,
    fitter = "glm",
    data = dd,
    values = list(X = seq(-1, 1, 0.1)),
    references = 0,
    contrasts = c("ratio", "difference")
  )
  expect_output(print(x))
})

test_that("standardize generates output with non-survival data (bootstrap)", {
  prob_predict.coxph <- function(...) 1 - riskRegression::predictRisk(...)
  require(survival)
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
    arguments = list(
      formula = Surv(U, D) ~ X + Z + X * Z,
      method = "breslow",
      x = TRUE,
      y = TRUE
    ),
    fitter = "coxph",
    data = dd,
    times = 1:5,
    predict_fun = prob_predict.coxph,
    values = list(X = seq(-1, 1, 0.1)),
    B = 10,
    references = 0,
    contrasts = "difference"
  )
  expect_output(print(x))
})

test_that("standardize generates output with non-survival data (no bootstrap)", {
  prob_predict.coxph <- function(...) 1 - riskRegression::predictRisk(...)
  require(survival)
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
    arguments = list(
      formula = Surv(U, D) ~ X + Z + X * Z,
      method = "breslow",
      x = TRUE,
      y = TRUE
    ),
    fitter = "coxph",
    data = dd,
    times = 1:5,
    predict_fun = prob_predict.coxph,
    values = list(X = seq(-1, 1, 0.1)),
    references = 0,
    contrasts = "difference"
  )
  expect_output(print(x))
})

test_that("standardize generates output with non-survival data (single time point)", {
  prob_predict.coxph <- function(...) 1 - riskRegression::predictRisk(...)
  require(survival)
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
    arguments = list(
      formula = Surv(U, D) ~ X + Z + X * Z,
      method = "breslow",
      x = TRUE,
      y = TRUE
    ),
    fitter = "coxph",
    data = dd,
    times = 3,
    predict_fun = prob_predict.coxph,
    values = list(X = seq(-1, 1, 0.1)),
    references = 0,
    contrasts = "difference"
  )
  expect_output(print(x))
})
