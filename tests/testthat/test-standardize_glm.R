test_that("bad names gives error", {
  set.seed(3)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  Y <- rnorm(n, mean = X + Z + 0.1 * X^2)
  dd <- data.frame(Z, X, Y)
  expect_error(standardize_glm(formula = Y ~ X + Z + I(X^2), data = dd, values = list(X = 1, blue = 2)))
})

test_that("non-binary exposure gives error when using the doubly robust estimator", {
  set.seed(4)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  Y <- rnorm(n, mean = X + Z + 0.1 * X^2)
  dd <- data.frame(Z, X, Y)
  expect_error(standardize_glm_dr(formula_outcome = Y ~ X + Z + I(X^2), formula_exposure = X ~ Z, fitbas, data = dd, values = list(X = c(-1, 0, 1))))
})

test_that("multiple exposure gives error when using the doubly robust estimator", {
  set.seed(5)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  X2 <- rnorm(n, mean = 0.5 * Z)
  Y <- rnorm(n, mean = X + Z + 0.1 * X^2)
  dd <- data.frame(Z, X, X2, Y)
  expect_error(standardize_glm_dr(formula_outcome = Y ~ X + X2 + Z + I(X^2), formula_exposure = X ~ Z, fitbas, data = dd, values = list(X = c(0, 1), X2 = c(0, 1))))
})

test_that("check that it fails with p_population unspecified", {
  singapore <- AF::singapore
  mi <- singapore$Age
  m <- mean(mi)
  s <- sd(mi)
  d <- 5
  expect_error(standardize_glm(
    formula = Oesophagealcancer ~ (Everhotbev + Age + Dial + Samsu + Cigs)^2,
    family = binomial, data = singapore,
    values = list(Everhotbev = 0:1), clusterid = "Set", case_control = TRUE,
    matched_density_cases = function(x) dnorm(x, m, s),
    matched_density_controls = function(x) dnorm(x, m - d, s), matching_variable = mi
  ))
})

test_that("non-valid transformation does not work", {
  set.seed(2)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  Y <- rnorm(n, mean = X + Z + 0.1 * X^2)
  dd <- data.frame(Z, X, Y)
  expect_error(standardize_glm(
    formula = Y ~ X * Z,
    family = "gaussian",
    data = dd,
    values = list(X = seq(-1, 1, 0.1)),
    transform = "log"
  ))
})

test_that("warning occurs when reference has been specified but not contrast", {
  data <- AF::clslowbwt
  expect_warning(standardize_glm_dr(
    formula_outcome = bwt ~ smoker * (race + age + lwt) + I(age^2) + I(lwt^2),
    formula_exposure = smoker ~ race * age * lwt + I(age^2) + I(lwt^2),
    family_outcome = gaussian,
    family_exposure = binomial,
    data = data,
    values = list(smoker = c(0, 1)), reference = 0
  ))
})

test_that("check estimates and standard errors standardize_glm (simple estimator)", {
  set.seed(6)
  n <- 100
  Z <- rnorm(n)
  X <- rnorm(n, mean = Z)
  Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
  dd <- data.frame(Z, X, Y)
  x <- standardize_glm(formula = Y ~ X * Z, family = "binomial", data = dd, values = list(X = 0:1))
  expect_equal(x$res_contrast[[1]]$estimates$estimates, c(0.519063874450474, 0.390531102199254), tolerance = 1e-5)
  expect_equal(x$res_contrast[[1]]$estimates$se, c(0.0614996028608024, 0.0881636202817664), tolerance = 1e-5)
})

test_that("check estimates and standard errors standardize_glm (case-control estimator)", {
  singapore <- AF::singapore
  mi <- singapore$Age
  m <- mean(mi)
  s <- sd(mi)
  d <- 5
  x <- suppressWarnings(standardize_glm(
    formula = Oesophagealcancer ~ (Everhotbev + Age + Dial + Samsu + Cigs)^2,
    family = binomial, data = singapore,
    values = list(Everhotbev = 0:1), clusterid = "Set", case_control = TRUE,
    matched_density_cases = function(x) dnorm(x, m, s),
    matched_density_controls = function(x) dnorm(x, m - d, s), matching_variable = mi,
    p_population = 19.3 / 100000
  ))
  expect_equal(x$res_contrast[[1]]$estimates$estimates, c(0.000127787218723152, 0.00057018810096904), tolerance = 1e-5)
  expect_equal(x$res_contrast[[1]]$estimates$se, c(1.68549810412155e-05, 0.000218937483527179), tolerance = 1e-5)
})

test_that("check estimates and standard errors standardize_glm (dr estimator)", {
  data <- AF::clslowbwt
  x <- standardize_glm_dr(
    formula_outcome = bwt ~ smoker * (race + age + lwt) + I(age^2) + I(lwt^2),
    formula_exposure = smoker ~ race * age * lwt + I(age^2) + I(lwt^2),
    family_outcome = gaussian,
    family_exposure = binomial,
    data = data,
    values = list(smoker = c(0, 1)), contrasts = "difference", reference = 0
  )
  expect_equal(x$res_contrast[[2]]$est_table$Estimate[2], -223.6736, tolerance = 1e-5)
  expect_equal(x$res_contrast[[2]]$est_table$`lower.0.95`[2], -424.6136, tolerance = 1e-5)
  expect_equal(x$res_contrast[[2]]$est_table$`upper.0.95`[2], -22.7335, tolerance = 1e-5)
})
