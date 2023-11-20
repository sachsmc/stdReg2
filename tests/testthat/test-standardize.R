test_that("standardize with non-survival data (no bootstrap) gives the same as standardize_glm", {
  set.seed(6)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
  dd <- data.frame(Z, X, Y)
  prob_predict.glm <- function(...) predict.glm(..., type = "response")
  y <- standardize_glm(
    formula = Y ~ X * Z,
    family = "binomial",
    data = dd,
    values = list(X = seq(-1, 1, 0.1)),
    reference = 0,
    contrasts = c("ratio", "difference")
  )
  x <- standardize(
    arguments = list(
      formula = Y ~ X * Z,
      family = "binomial"
    ),
    predict_fun = prob_predict.glm,
    fitter = "glm",
    data = dd,
    values = list(X = seq(-1, 1, 0.1)),
    reference = 0,
    contrasts = c("ratio", "difference")
  )
  expect_equal(x$res$estimates[, 2], y$res$estimates[, 3])
})

test_that("standardize with survival data (no bootstrap) gives the same as standardize_coxph", {
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
  time <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
  C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
  U <- pmin(time, C) # time at risk
  D <- as.numeric(time < C) # event indicator
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
    reference = 0,
    contrasts = "difference"
  )
  fit.std <- standardize_coxph(
    formula = Surv(U, D) ~ X + Z + X * Z,
    data = dd,
    values = list(X = seq(-1, 1, 0.1)),
    times = 1:5,
    reference = 0,
    contrasts = "difference"
  )

  expect_equal(fit.std$res$est, unname(t(x$res$estimates[, -1])))

})

test_that("standardize generates output with survival data (single time point)", {
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
  time <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
  C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
  U <- pmin(time, C) # time at risk
  D <- as.numeric(time < C) # event indicator
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
    reference = 0,
    contrasts = "difference"
  )
  expect_output(print(x))
})
