#' @title Get standardized estimates using the g-formula with a custom model
#' @param arguments
#' The arguments to be used in the fitter function as a \code{list}.
#' @param fitter The function to call to fit the data.
#' @param predict_fun The function used to predict the means/probabilities
#' for a new data set on the response level. For survival data,
#' this should be a matrix where each column is the time, and each
#' row the data.
#' @param times For use with survival data. Set to \code{NULL} otherwise.
#' @param data
#' The data used.
#' @param values
#' A named list or data.frame specifying the variables and values
#' at which marginal means of the outcome will be estimated.
#' @returns
#' An object of class \code{std_helper}.
#' This is basically a list with components estimates and fit for the outcome model.
#' @details
#' Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome, the exposure, and a
#' vector of covariates, respectively.
#' \code{standardize} uses a
#' model to estimate the standardized
#' mean \eqn{\theta(x)=E\{E(Y|X=x,Z)\}},
#' where \eqn{x} is a specific value of \eqn{X},
#' and the outer expectation is over the marginal distribution of \eqn{Z}.
#' With survival data, \eqn{Y=I(T > t)},
#' and a vector of different time points \code{times} (\eqn{t}) can be given,
#' where \eqn{T} is the uncensored survival time.
#' @examples
#'
#' set.seed(6)
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
#' dd <- data.frame(Z, X, Y)
#' prob_predict.glm <- function(...) predict.glm(..., type = "response")
#'
#' x <- standardize(
#'   arguments = list(
#'     formula = Y ~ X * Z,
#'     family = "binomial"
#'   ),
#'   predict_fun = prob_predict.glm,
#'   fitter = "glm",
#'   data = dd,
#'   values = list(X = seq(-1, 1, 0.1))
#' )$estimates
#'
#' x
#'
#' prob_predict.coxph <- function(...) 1 - riskRegression::predictRisk(...)
#' require(survival)
#' set.seed(68)
#' n <- 500
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
#' C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
#' U <- pmin(T, C) # time at risk
#' D <- as.numeric(T < C) # event indicator
#' dd <- data.frame(Z, X, U, D)
#' x <- standardize(
#'   arguments = list(
#'     formula = Surv(U, D) ~ X + Z + X * Z,
#'     method = "breslow",
#'     x = TRUE,
#'     y = TRUE
#'   ),
#'   fitter = "coxph",
#'   data = dd,
#'   times = 1:5,
#'   predict_fun = prob_predict.coxph,
#'   values = list(X = seq(-1, 1, 0.1))
#' )$estimates
#' x
#'
#' @export standardize
standardize <- function(arguments,
                        fitter,
                        times = NULL,
                        predict_fun,
                        data,
                        values) {
  ## Preparation and various checks
  n <- nrow(data)

  if (!inherits(values, c("data.frame", "list"))) {
    stop("values is not an object of class list or data.frame")
  }

  ## Check that the names of values appear in the data
  check_values_data(values, data)

  ## Set various relevant variables
  if (!is.data.frame(values)) {
    valuesout <- expand.grid(values)
  } else {
    valuesout <- values
  }
  exposure_names <- colnames(valuesout)
  exposure <- data[, exposure_names]

  fit_outcome <- fit_helper(arguments, fitter, data)
  if (is.null(times)) {
    estimates <- rep(NA, nrow(valuesout))
  } else {
    estimates <- matrix(NA, ncol = length(times), nrow = nrow(valuesout))
    colnames(estimates) <- paste("t =", times)
  }

  for (i in seq_len(nrow(valuesout))) {
    data_x <- do.call("transform", c(
      list(data),
      valuesout[i, , drop = FALSE]
    ))

    ## Save the predictions for data_x
    if (is.null(times)) {
      estimates[i] <- mean(predict_fun(object = fit_outcome, newdata = data_x))
    } else {
      estimates[i, ] <- colMeans(predict_fun(object = fit_outcome, newdata = data_x, times = times))
    }
  }
  valuesout <- cbind(valuesout, estimates)
  res <- structure(
    list(
      estimates = valuesout,
      fit_outcome = fit_outcome
    ),
    class = "std_helper"
  )
  res
}

#' @title Get the causal mean restricted survival time using the g-formula with a custom model
#' @param time_grid The time grid to be used for numerical integration.
#' @inherit standardize
#' @details
#' Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome, the exposure, and a
#' vector of covariates, respectively.
#' \code{standardize} uses a
#' model to estimate the standardized
#' mean \eqn{\theta(x)=\int_0^L E\{E(I(T > u)|X=x,Z)\}} du,
#' where \eqn{x} is a specific value of \eqn{X},
#' and the outer expectation is over the marginal distribution of \eqn{Z}, where \eqn{T}
#' is the uncensored survival time. Here \eqn{L} is end-of-study specified as the last
#' time point in \code{time_grid}. The integral is approximated by a \code{time_grid}.
#' @examples
#'
#' library(survival)
#' prob_predict.coxph <- function(...) 1 - riskRegression::predictRisk(...)
#' set.seed(6)
#' n <- 500
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
#' C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
#' U <- pmin(T, C) # time at risk
#' D <- as.numeric(T < C) # event indicator
#' dd <- data.frame(Z, X, U, D)
#' x <- standardize_restricted_mean(
#'   arguments = list(
#'     formula = Surv(U, D) ~ X + Z + X * Z,
#'     method = "breslow",
#'     x = TRUE,
#'     y = TRUE
#'   ),
#'   fitter = "coxph",
#'   data = dd,
#'   time_grid = seq(0, 5, 0.001)[-1],
#'   predict_fun = prob_predict.coxph,
#'   values = list(X = 0.5)
#' )
#' @export standardize_restricted_mean
standardize_restricted_mean <- function(arguments,
                                        fitter,
                                        time_grid,
                                        predict_fun,
                                        data,
                                        values) {
  ## Preparation and various checks
  n <- nrow(data)

  if (!inherits(values, c("data.frame", "list"))) {
    stop("values is not an object of class list or data.frame")
  }

  ## Check that the names of values appear in the data
  check_values_data(values, data)

  ## Set various relevant variables
  if (!is.data.frame(values)) {
    valuesout <- expand.grid(values)
  } else {
    valuesout <- values
  }
  exposure_names <- colnames(valuesout)
  exposure <- data[, exposure_names]

  fit_outcome <- fit_helper(arguments, fitter, data)
  estimates <- rep(NA, nrow(valuesout))

  for (i in seq_len(nrow(valuesout))) {
    data_x <- do.call("transform", c(
      list(data),
      valuesout[i, , drop = FALSE]
    ))
    temp <- predict_fun(object = fit_outcome, newdata = data_x, times = time_grid)
    trapezoidal_rule <- function(time_grid, y_val) {
      sum((y_val[1:(length(time_grid) - 1)] + y_val[2:(length(time_grid))]) / 2 * diff(time_grid))
    }
    estimates[i] <- mean(apply(temp, 1, function(y_val) trapezoidal_rule(time_grid, y_val)))
  }
  valuesout <- cbind(valuesout, estimates)
  res <- structure(
    list(
      estimates = valuesout,
      fit_outcome = fit_outcome
    ),
    class = "std_helper"
  )
  res
}

#' @rdname print
#' @export print.std_helper
#' @export
print.std_helper <- function(x, ...) {
  print(x$estimates)
}


fit_helper <- function(args, fitter, data) {
  ## try fitting a glm model
  args[["data"]] <- data
  fit <- tryCatch(
    {
      do.call(fitter, args)
    },
    error = function(cond) {
      return(cond)
    }
  )
  if (inherits(fit, "simpleError")) {
    stop("fitter failed with error: ", fit[["message"]])
  } else {
    fit
  }
}
