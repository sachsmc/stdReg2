#' @title Regression standardization in Cox proportional hazards models
#'
#' @description \code{standardize_coxph} performs regression standardization in Cox proportional
#' hazards models, at specified values of the exposure, over the sample
#' covariate distribution. Let \eqn{T}, \eqn{X}, and \eqn{Z} be the survival
#' outcome, the exposure, and a vector of covariates, respectively.
#' \code{standardize_coxph} uses a fitted Cox proportional hazards model to estimate the
#' standardized survival function \eqn{\theta(t,x)=E\{S(t|X=x,Z)\}}, where
#' \eqn{t} is a specific value of \eqn{T}, \eqn{x} is a specific value of
#' \eqn{X}, and the expectation is over the marginal distribution of \eqn{Z}.
#'
#' @details \code{standardize_coxph} assumes that a Cox proportional hazards model
#' \deqn{\lambda(t|X,Z)=\lambda_0(t)exp\{h(X,Z;\beta)\}} has been fitted.
#' Breslow's estimator of the cumulative baseline hazard
#' \eqn{\Lambda_0(t)=\int_0^t\lambda_0(u)du} is used together with the partial
#' likelihood estimate of \eqn{\beta} to obtain estimates of the survival
#' function \eqn{S(t|X=x,Z)}:
#' \deqn{\hat{S}(t|X=x,Z)=exp[-\hat{\Lambda}_0(t)exp\{h(X=x,Z;\hat{\beta})\}].}
#' For each \eqn{t} in the \code{t} argument and for each \eqn{x} in the
#' \code{x} argument, these estimates are averaged across all subjects (i.e.
#' all observed values of \eqn{Z}) to produce estimates
#' \deqn{\hat{\theta}(t,x)=\sum_{i=1}^n \hat{S}(t|X=x,Z_i)/n,} where \eqn{Z_i}
#' is the value of \eqn{Z} for subject \eqn{i}, \eqn{i=1,...,n}.  The variance
#' for \eqn{\hat{\theta}(t,x)} is obtained by the sandwich formula.
#' @return An object of class \code{std_surv}.
#' This is basically a list with components estimates and covariance estimates in \code{res}
#' Results for transformations, contrasts, references are stored in \code{res_contrasts}.
#'  The output contains estimates for contrasts and confidence intervals for all combinations of transforms, references
#' @inherit standardize_glm
#' @param times A vector containing the specific values of \eqn{T} at
#' which to estimate the standardized survival function.
#' @note Standardized survival functions are sometimes referred to as (direct)
#' adjusted survival functions in the literature.
#'
#' \code{standardize_coxph/standardize_parfrailty} does not currently handle time-varying exposures or
#' covariates.
#'
#' \code{standardize_coxph/standardize_parfrailty} internally loops over all values in the \code{t} argument.
#' Therefore, the function will usually be considerably faster if
#' \code{length(t)} is small.
#'
#' The variance calculation performed by \code{standardize_coxph} does not condition on
#' the observed covariates \eqn{\bar{Z}=(Z_1,...,Z_n)}. To see how this
#' matters, note that
#' \deqn{var\{\hat{\theta}(t,x)\}=E[var\{\hat{\theta}(t,x)|\bar{Z}\}]+var[E\{\hat{\theta}(t,x)|\bar{Z}\}].}
#' The usual parameter \eqn{\beta} in a Cox proportional hazards model does not
#' depend on \eqn{\bar{Z}}. Thus, \eqn{E(\hat{\beta}|\bar{Z})} is independent
#' of \eqn{\bar{Z}} as well (since \eqn{E(\hat{\beta}|\bar{Z})=\beta}), so that
#' the term \eqn{var[E\{\hat{\beta}|\bar{Z}\}]} in the corresponding variance
#' decomposition for \eqn{var(\hat{\beta})} becomes equal to 0. However,
#' \eqn{\theta(t,x)} depends on \eqn{\bar{Z}} through the average over the
#' sample distribution for \eqn{Z}, and thus the term
#' \eqn{var[E\{\hat{\theta}(t,x)|\bar{Z}\}]} is not 0, unless one conditions on
#' \eqn{\bar{Z}}. The variance calculation by Gail and Byar (1986) ignores this
#' term, and thus effectively conditions on \eqn{\bar{Z}}.
#' @author Arvid Sjolander
#' @references
#'
#' Chang I.M., Gelman G., Pagano M. (1982). Corrected group prognostic curves
#' and summary statistics. \emph{Journal of Chronic Diseases} \bold{35},
#' 669-674.
#'
#' Gail M.H. and Byar D.P. (1986). Variance calculations for direct adjusted
#' survival curves, with applications to testing for no treatment effect.
#' \emph{Biometrical Journal} \bold{28}(5), 587-599.
#'
#' Makuch R.W. (1982). Adjusted survival curve estimation using covariates.
#' \emph{Journal of Chronic Diseases} \bold{35}, 437-443.
#'
#' Sjolander A. (2016). Regression standardization with the R-package stdReg.
#' \emph{European Journal of Epidemiology} \bold{31}(6), 563-574.
#'
#' Sjolander A. (2016). Estimation of causal effect measures with the R-package
#' stdReg. \emph{European Journal of Epidemiology} \bold{33}(9), 847-858.
#' @examples
#'
#'
#' require(survival)
#' set.seed(7)
#' n <- 300
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
#' C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
#' U <- pmin(T, C) # time at risk
#' D <- as.numeric(T < C) # event indicator
#' dd <- data.frame(Z, X, U, D)
#' fit.std <- standardize_coxph(
#'   formula = Surv(U, D) ~ X + Z + X * Z,
#'   data = dd,
#'   values = list(X = seq(-1, 1, 0.5)),
#'   times = 1:5
#' )
#' print(fit.std)
#' plot(fit.std)
#'
#' @export standardize_coxph
standardize_coxph <- function(formula,
                              data,
                              values,
                              times,
                              clusterid,
                              ci_level = 0.95,
                              ci_type = "plain",
                              contrasts = NULL,
                              family = "gaussian",
                              references = NULL,
                              transforms = NULL) {
  call <- match.call()

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

  fit <- tryCatch(
    {
      survival::coxph(formula = formula, data = data, method = "breslow", x = TRUE, y = TRUE)
    },
    error = function(cond) {
      return(cond)
    }
  )
  if (inherits(fit, "simpleError")) {
    stop("fitter failed with error: ", fit[["message"]])
  }

  #---PREPARATION---
  specials <- pmatch(c("strata(", "cluster(", "tt("), attr(
    terms(fit$formula),
    "variables"
  ))
  if (any(!is.na(specials))) {
    stop("No special terms are allowed in the formula")
  }

  npar <- length(fit$coef)
  fit.detail <- coxph.detail(object = fit)

  # Delete rows that did not contribute to the model fit,
  # e.g. missing data or not in subset for fit.
  # Need to have object=fit in model.matrix, since neither object=formula nor
  # object=terms(fit) will remove rows not in subset.
  m <- model.matrix(object = fit)
  data <- data[match(rownames(m), rownames(data)), ]
  n <- nrow(data)

  input <- as.list(environment())

  if (is.null(fit$weights)) {
    weights <- rep(1, nrow(data))
  } else {
    weights <- fit$weights
  }

  # Can write code more generally with
  # if(missing(clusters)) clusters <- 1:nrow(data)
  # but a problem when constructing meat in sandwich formula:
  # must always aggregate, which is a bit slow, even though much faster
  # when using data.table than the aggregate function.
  if (!missing(clusterid)) {
    ncluster <- length(unique(data[, clusterid]))
  }

  nX <- nrow(valuesout)

  # Assign value to times if missing.
  if (missing(times)) {
    stop("You have to specify the times (t) at which to estimate the standardized survival function")
    # t <- fit.detail$time
  }
  input$times <- times
  nt <- length(times)

  if (sum(fit.detail$time <= min(times)) == 0) {
    stop("No events before first value in times", call. = FALSE)
  }

  est <- matrix(nrow = nt, ncol = nX)
  vcov <- vector(mode = "list", length = nt)
  H <- Hfun(fit = fit, data = data, fit.detail = fit.detail)
  sandwich.fit <- sandwich(
    fit = fit, data = data, weights = weights, t = times,
    fit.detail = fit.detail
  )

  #---LOOP OVER nt

  for (j in 1:nt) {
    if (times[j] == 0) {
      est[j, ] <- 1
      vcov[[j]] <- matrix(0, nrow = nX, ncol = nX)
    } else {
      #---ESTIMATES OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---

      si <- matrix(nrow = n, ncol = nX)
      PredX <- matrix(nrow = n, ncol = nX)
      tempmat <- matrix(nrow = nX, ncol = npar)
      for (i in 1:nX) {
        data.x <- do.call("transform", c(
          list(data),
          valuesout[i, , drop = FALSE]
        ))
        predX <- predict(object = fit, newdata = data.x, type = "risk")
        si[, i] <- exp(-H(times[j]) * predX)
        PredX[, i] <- predX
        # Need terms(fit) here. If formula contains splines,
        # then fit or formula will not work when no variation in the exposure,
        # since model.matrix need to retrieve Boundary.knots from terms(fit).
        # Also need to center, since everything else is centered.
        # Note: need to center at factual X, not counterfactual x,
        # since baseline is computed at mean of factual X.
        m.x <- model.matrix(object = terms(fit), data = data.x)[, -1, drop = FALSE]
        m <- model.matrix(object = terms(fit), data = data)[, -1, drop = FALSE]
        m <- matrix(colMeans(m), nrow = nrow(m), ncol = ncol(m), byrow = TRUE)
        m.x <- m.x - m
        tempmat[i, ] <- colMeans(m.x * predX * si[, i] * weights)
      }
      est[j, ] <- colSums(weights * si, na.rm = TRUE) /
        sum(weights)

      #---VARIANCE OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x, ---

      sres <- weights * (si - matrix(rep(est[j, ], each = n), nrow = n, ncol = nX))
      ores <- sandwich.fit$U[, c(1:npar, npar + j)]
      res <- cbind(sres, ores)
      if (!missing(clusterid)) {
        res <- aggr(x = res, clusters = data[, clusterid])
      }
      J <- var(res, na.rm = TRUE)
      SI <- cbind(
        -diag(nX) * mean(weights), -tempmat * H(times[j]),
        -colMeans(PredX * si * weights)
      )
      # This is why the user cannot use term cluster; then -solve(vcov(object=fit))/n
      # will not be the bread in the sandwich.
      oI <- cbind(
        matrix(0, nrow = npar + 1, ncol = nX),
        sandwich.fit$I[c(1:npar, npar + j), c(1:npar, npar + j)]
      )
      I <- rbind(SI, oI)

      if (missing(clusterid)) {
        V <- (solve(I) %*% J %*% t(solve(I)) / n)[1:nX, 1:nX]
      } else {
        V <- (solve(I) %*% J %*% t(solve(I)) * ncluster / n^2)[1:nX, 1:nX]
      }
      vcov[[j]] <- V
    }
  }

  out <- list(call = call, input = input, est = est, vcov = vcov)
  class(out) <- "std_coxph"

  #---OUTPUT---
  format_result_standardize(
    out,
    contrasts,
    references,
    transforms,
    ci_type,
    ci_level,
    "std_surv",
    "summary_std_coxph"
  )
}

summary_std_coxph <- function(object, times, ci_type = "plain", ci_level = 0.95,
                              transform = NULL, contrast = NULL, reference = NULL, ...) {
  if (dim(object$input$valuesout)[2] > 1) {
    stop("multiple exposures not currently suported with standardize_coxph/standardize_parfrailty")
  }
  est.all <- object$est
  V.all <- object$vcov
  nX <- nrow(object$input$valuesout)
  if (missing(times)) {
    times <- object$input$times
  }
  nt <- length(times)

  est.table <- vector(mode = "list", length = nt)
  for (j in 1:nt) {
    if (min(abs(times[j] - object$input$times)) > sqrt(.Machine$double.eps)) {
      stop("The standardized survival function is not estimated at times",
        call. = FALSE
      )
    } else {
      k <- which.min(abs(times[j] - object$input$times))
    }

    est <- est.all[k, ]
    V <- as.matrix(V.all[[k]])

    if (!is.null(transform)) {
      if (transform == "log") {
        dtransform.dm <- diag(1 / est, nrow = nX, ncol = nX)
        est <- log(est)
      }
      if (transform == "logit") {
        dtransform.dm <- diag(1 / (est * (1 - est)), nrow = nX, ncol = nX)
        est <- logit(est)
      }
      if (transform == "odds") {
        dtransform.dm <- diag(1 / (1 - est)^2, nrow = nX, ncol = nX)
        est <- odds(est)
      }
      V <- t(dtransform.dm) %*% V %*% dtransform.dm
    }

    if (!is.null(contrast)) {
      if (is.null(reference)) {
        stop("When specifying contrast, reference must be specified as well")
      }
      referencepos <- match(reference, object$input$valuesout[, 1])
      if (is.na(referencepos)) {
        stop("reference must be a value in x")
      }
      if (contrast == "difference") {
        dcontrast.dtransform <- diag(nX)
        dcontrast.dtransform[referencepos, ] <- -1
        dcontrast.dtransform[referencepos, referencepos] <- 0
        est <- est - est[referencepos]
      }
      if (contrast == "ratio") {
        dcontrast.dtransform <- diag(1 / est[referencepos], nrow = nX, ncol = nX)
        dcontrast.dtransform[referencepos, ] <- -est / est[referencepos]^2
        dcontrast.dtransform[referencepos, referencepos] <- 1
        est <- est / est[referencepos]
      }
      V <- t(dcontrast.dtransform) %*% V %*% dcontrast.dtransform
      V[referencepos, ] <- 0
      V[, referencepos] <- 0
    }

    var <- diag(V)
    se <- sqrt(var)
    conf.int <- CI(est = est, var = var, ci_type = ci_type, ci_level = ci_level)

    temp <- as.matrix(cbind(est, se, conf.int), nrow = length(est), ncol = 4)
    dimnames(temp) <- list(
      object$input$x,
      c(
        "Estimate", "Std. Error", paste("lower", ci_level),
        paste("upper", ci_level)
      )
    )
    est.table[[j]] <- temp
  }
  if (is.factor(reference)) {
    reference <- as.character(reference)
  }
  out <- c(
    object,
    list(
      est.table = est.table, tsum = times, transform = transform, contrast = contrast,
      reference = reference
    )
  )
  return(out)
}

#' @rdname print
#' @export print.std_surv
#' @export
print.std_surv <- function(x, ...) {
  for (v in x$res_contrast) {
    print_summary_std_coxph(summary_std_coxph(v))
  }
}

print_summary_std_coxph <- function(x, ...) {
  nt <- length(x$tsum)
  for (j in 1:nt) {
    cat("\nFormula: ")
    print(x$input$fit$formula, showEnv = FALSE)
    cat("Exposure: ", x$input$exposure_names, "\n")

    if (!is.null(x$transform)) {
      cat("Transform: ", x$transform, "\n")
    }
    if (!is.null(x$contrast)) {
      cat("Reference level: ", x$input$X, "=", x$reference, "\n")
      cat("Contrast: ", x$contrast, "\n")
    }
    cat("Survival functions evaluated at t =", x$tsum[j], "\n")
    cat("\n")
    print(x$est.table[[j]], digits = 3)
    cat("\n")
  }
}

#' @param legendpos position of the legend; see help for \code{legend}.
#' @rdname plot
#' @export plot.std_surv
#' @export
plot.std_surv <- function(x, plot_ci = TRUE, ci_type = "plain", ci_level = 0.95,
                          transform = NULL, contrast = NULL, reference = NULL, legendpos = "bottomleft", summary_fun = "summary_std_coxph", ...) {
  object <- x$res
  if (ncol(object$input$valuesout) != 1) {
    stop("multiple exposures")
  } else {
    x <- object$input$valuesout[, 1]
  }

  dots <- list(...)

  xlab <- "t"

  if (is.factor(reference)) {
    reference <- as.character(reference)
  }

  if (is.null(contrast)) {
    if (is.null(transform)) {
      ylab <- expression(S(t))
    } else {
      if (transform == "log") {
        ylab <- expression(paste(log, "{", S(t), "}", sep = ""))
      }
      if (transform == "logit") {
        ylab <- expression(paste(logit, "{", S(t), "}", sep = ""))
      }
      if (transform == "odds") {
        ylab <- expression(paste(S(t), "/{", 1 - S(t), "}", sep = ""))
      }
    }
  } else {
    if (contrast == "difference") {
      if (is.null(transform)) {
        ylab <- c(bquote(paste(S(t), "-", S[.(reference)](t))), expression())
      } else {
        if (transform == "log") {
          ylab <- c(bquote(paste(log, "{", S(t), "}-", log, "{",
            S[.(reference)](t), "}",
            sep = ""
          )), expression())
        }
        if (transform == "logit") {
          ylab <- c(bquote(paste(logit, "{", S(t), "}-", logit,
            "{", S[.(reference)](t), "}",
            sep = ""
          )), expression())
        }
        if (transform == "odds") {
          ylab <- c(bquote(paste(S(t), "/{", 1 - S(t), "}-",
            S[.(reference)](t), "/{", 1 - S[.(reference)](t),
            "}",
            sep = ""
          )), expression())
        }
      }
    }
    if (contrast == "ratio") {
      if (is.null(transform)) {
        ylab <- c(
          bquote(paste(S(t), " / ", S[.(reference)](t), sep = "")),
          expression()
        )
      } else {
        if (transform == "log") {
          ylab <- c(bquote(paste(log, "{", S(t), "} / ", log,
            "{", S[.(reference)](t), "}",
            sep = ""
          )), expression())
        }
        if (transform == "logit") {
          ylab <- c(bquote(paste(logit, "{", S(t), "} / ", logit,
            "{", S[.(reference)](t), "}",
            sep = ""
          )), expression())
        }
        if (transform == "odds") {
          ylab <- c(bquote(paste("[", S(t), "/{", 1 - S(t), "}] / [",
            S[.(reference)](t), "/{", 1 - S[.(reference)](t),
            "}]",
            sep = ""
          )), expression())
        }
      }
    }
  }

  times <- object$input$times
  nt <- length(times)
  nX <- length(x)

  sum.obj <- do.call(summary_fun, list(
    object = object, ci_type = ci_type, ci_level = ci_level,
    transform = transform, contrast = contrast, reference = reference
  ))

  temp <- Reduce(f = rbind, x = sum.obj$est.table)
  est <- matrix(temp[, 1], nrow = nt, ncol = nX, byrow = TRUE)
  lower <- matrix(temp[, 3], nrow = nt, ncol = nX, byrow = TRUE)
  upper <- matrix(temp[, 4], nrow = nt, ncol = nX, byrow = TRUE)
  if (plot_ci) {
    ylim <- c(min(lower), max(upper))
  } else {
    ylim <- c(min(est), max(est))
  }
  args <- list(
    x = times, y = rep(0, length(times)), xlab = xlab, ylab = ylab,
    ylim = ylim, type = "n"
  )
  args[names(dots)] <- dots
  do.call("plot", args = args)
  legend <- NULL
  for (i in seq_len(nX)) {
    lines(times, est[, i], col = i)
    if (plot_ci) {
      lines(times, upper[, i], lty = "dashed", col = i)
      lines(times, lower[, i], lty = "dashed", col = i)
    }
    temp <- as.character(x[i])
  }
  legend <- c(legend, paste(object$input$exposure_names, "=", object$input$valuesout[, 1]))
  legend(
    x = legendpos, legend = legend, lty = rep(1, length(x)), col = seq_len(length(x)),
    bty = "n"
  )
}
