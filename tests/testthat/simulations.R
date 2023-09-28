# library(stdReg)
# library(emmeans)
# library(marginaleffects)
#
#
#
# dat <- mtcars
# dat$carb <- factor(dat$carb)
# dat$cyl <- factor(dat$cyl)
# dat$am <- as.logical(dat$am)
# mod <- glm(mpg ~ carb + cyl + am, dat, family = "gaussian")
#
# marginal_means(
#   mod,
#   variables = "cyl"
# )
#
# fit.std <- stdGlm(fit = mod, data = dat, X = "cyl", x = c("4", "6", "8"))
# print(summary(fit.std))
#
# emmeans(mod, "cyl")
# mean(predict(mod, newdata = subset(dat, cyl == "4")))
#
# chke <- ref_grid(mod)@grid
# chke$mean <- predict(mod, newdata = ref_grid(mod)@grid)
# by(chke$mean, chke$cyl, mean)
#
# stdgrid <- ref_grid(mod, counterfactuals = "cyl")
# summary(stdgrid)
#
#
# ## true values
#
# theta <- -integrate(\(xx) plogis(1 + xx + 1 * xx) * dnorm(xx, sd = 4), -Inf, Inf)$value +
#   integrate(\(xx) plogis(xx) * dnorm(xx, sd = 4), -Inf, Inf)$value
#
# simswag <- replicate(2000, {
#   n <- 500
#   Z <- rnorm(n, sd = 4)
#   X <- rbinom(n, 1, prob = plogis(Z))
#   Y <- rbinom(n, 1, prob = plogis(X + Z + X * Z))
#
#   dd <- data.frame(Z, X, Y)
#   fit <- glm(formula = Y ~ X + Z + X * Z, family = "binomial", data = dd)
#   fit.std <- stdGlm(fit = fit, data = dd, X = "X", x = 0:1)
#   cis <- list(
#     summary(fit.std, contrast = "difference", reference = 1)$est.table[1, 3:4],
#     confint(pairs(ref_grid(fit, counterfactuals = "X"))) |> (\(x){
#       c(x$asymp.LCL, x$asymp.UCL)
#     })()
#   )
#
#   sapply(cis, \(x) {
#     theta >= x[1] & theta <= x[2]
#   }) |> setNames(c("stdreg", "emmeans"))
# })
#
# rowMeans(simswag)
#
#
# thetar <- integrate(\(xx) plogis(1 + xx + 1 * xx) * dnorm(xx, sd = 4), -Inf, Inf)$value /
#   integrate(\(xx) plogis(xx) * dnorm(xx, sd = 4), -Inf, Inf)$value
#
# simswag <- replicate(2000, {
#   n <- 500
#   Z <- rnorm(n, sd = 4)
#   X <- rbinom(n, 1, prob = plogis(Z))
#   Y <- rbinom(n, 1, prob = plogis(X + Z + X * Z))
#
#   dd <- data.frame(Z, X, Y)
#   fit <- glm(formula = Y ~ X + Z + X * Z, family = "binomial", data = dd)
#   fit.std <- stdGlm(fit = fit, data = dd, X = "X", x = 0:1)
#   cis <- list(
#     summary(fit.std, contrast = "ratio", reference = 0)$est.table[2, 3:4],
#     confint(ref_grid(fit, counterfactuals = "X") |>
#       regrid(transform = "log") |>
#       contrast("trt.vs.ctrl"), type = "response") |>
#       (\(x){
#         c(x$asymp.LCL, x$asymp.UCL)
#       })()
#   )
#
#   sapply(cis, \(x) {
#     thetar >= x[1] & thetar <= x[2]
#   }) |> setNames(c("stdreg", "emmeans"))
# })
#
#
# rowMeans(simswag)
