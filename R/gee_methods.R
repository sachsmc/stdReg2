#' @title Regression standardization in conditional generalized estimating equations
#'
#' @description \code{standardize_gee} performs regression standardization in linear and log-linear
#' fixed effects models, at specified values of the exposure, over the sample
#' covariate distribution. Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome,
#' the exposure, and a vector of covariates, respectively. It is assumed that
#' data are clustered with a cluster indicator \eqn{i}. \code{standardize_gee} uses
#' fitted fixed effects model, with cluster-specific intercept \eqn{a_i} (see
#' \code{details}), to estimate the standardized mean
#' \eqn{\theta(x)=E\{E(Y|i,X=x,Z)\}}, where \eqn{x} is a specific value of
#' \eqn{X}, and the outer expectation is over the marginal distribution of
#' \eqn{(a_i,Z)}.
#'
#' @details \code{standardize_gee} assumes that a fixed effects model
#' \deqn{\eta\{E(Y|i,X,Z)\}=a_i+h(X,Z;\beta)} has been fitted. The link
#' function \eqn{\eta} is assumed to be the identity link or the log link. The
#' conditional generalized estimating equation (CGEE) estimate of \eqn{\beta}
#' is used to obtain estimates of the cluster-specific means:
#' \deqn{\hat{a}_i=\sum_{j=1}^{n_i}r_{ij}/n_i,} where
#' \deqn{r_{ij}=Y_{ij}-h(X_{ij},Z_{ij};\hat{\beta})} if \eqn{\eta} is the
#' identity link, and \deqn{r_{ij}=Y_{ij}\exp\{-h(X_{ij},Z_{ij};\hat{\beta})\}}
#' if \eqn{\eta} is the log link, and \eqn{(X_{ij},Z_{ij})} is the value of
#' \eqn{(X,Z)} for subject \eqn{j} in cluster \eqn{i}, \eqn{j=1,...,n_i},
#' \eqn{i=1,...,n}. The CGEE estimate of \eqn{\beta} and the estimate of
#' \eqn{a_i} are used to estimate the mean \eqn{E(Y|i,X=x,Z)}:
#' \deqn{\hat{E}(Y|i,X=x,Z)=\eta^{-1}\{\hat{a}_i+h(X=x,Z;\hat{\beta})\}.} For
#' each \eqn{x} in the \code{x} argument, these estimates are averaged across
#' all subjects (i.e. all observed values of \eqn{Z} and all estimated values
#' of \eqn{a_i}) to produce estimates \deqn{\hat{\theta}(x)=\sum_{i=1}^n
#' \sum_{j=1}^{n_i} \hat{E}(Y|i,X=x,Z_i)/N,} where \eqn{N=\sum_{i=1}^n n_i}.
#' The variance for \eqn{\hat{\theta}(x)} is obtained by the sandwich formula.
#'
#' @param formula A formula to be used with \code{"gee"} in the \pkg{drgee} package.
#' @param link The link function to be used with \code{"gee"}.
#' @inherit standardize_glm
#' @note The variance calculation performed by \code{standardize_gee} does not condition
#' on the observed covariates \eqn{\bar{Z}=(Z_{11},...,Z_{nn_i})}. To see how
#' this matters, note that
#' \deqn{var\{\hat{\theta}(x)\}=E[var\{\hat{\theta}(x)|\bar{Z}\}]+var[E\{\hat{\theta}(x)|\bar{Z}\}].}
#' The usual parameter \eqn{\beta} in a generalized linear model does not
#' depend on \eqn{\bar{Z}}. Thus, \eqn{E(\hat{\beta}|\bar{Z})} is independent
#' of \eqn{\bar{Z}} as well (since \eqn{E(\hat{\beta}|\bar{Z})=\beta}), so that
#' the term \eqn{var[E\{\hat{\beta}|\bar{Z}\}]} in the corresponding variance
#' decomposition for \eqn{var(\hat{\beta})} becomes equal to 0. However,
#' \eqn{\theta(x)} depends on \eqn{\bar{Z}} through the average over the sample
#' distribution for \eqn{Z}, and thus the term
#' \eqn{var[E\{\hat{\theta}(x)|\bar{Z}\}]} is not 0, unless one conditions on
#' \eqn{\bar{Z}}.
#' @author Arvid Sjölander.
#' @references Goetgeluk S. and Vansteelandt S. (2008). Conditional generalized
#' estimating equations for the analysis of clustered and longitudinal data.
#' \emph{Biometrics} \bold{64}(3), 772-780.
#'
#' Martin R.S. (2017). Estimation of average marginal effects in multiplicative
#' unobserved effects panel models. \emph{Economics Letters} \bold{160}, 16-19.
#'
#' Sjölander A. (2019). Estimation of marginal causal effects in the presence
#' of confounding by cluster. \emph{Biostatistics} doi:
#' 10.1093/biostatistics/kxz054
#' @examples
#'
#' require(drgee)
#'
#' set.seed(4)
#' n <- 300
#' ni <- 2
#' id <- rep(1:n, each = ni)
#' ai <- rep(rnorm(n), each = ni)
#' Z <- rnorm(n * ni)
#' X <- rnorm(n * ni, mean = ai + Z)
#' Y <- rnorm(n * ni, mean = ai + X + Z + 0.1 * X^2)
#' dd <- data.frame(id, Z, X, Y)
#' fit.std <- standardize_gee(
#'   formula = Y ~ X + Z + I(X^2),
#'   link = "identity",
#'   data = dd,
#'   values = list(X = seq(-3, 3, 0.5)),
#'   clusterid = "id"
#' )
#' print(fit.std)
#' plot(fit.std)
#'
#' @export standardize_gee
standardize_gee <- function(formula, link = "identity", data, values, clusterid,
                            case_control = FALSE,
                            ci_level = 0.95,
                            ci_type = "plain",
                            contrasts = NULL,
                            family = "gaussian",
                            reference = NULL,
                            transforms = NULL) {
  if (!inherits(values, c("data.frame", "list"))) {
    stop("values is not an object of class list or data.frame")
  }

  n <- nrow(data)

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

  fit <- tryCatch(
    {
      drgee::gee(formula = formula, data = data, cond = TRUE, clusterid = clusterid, link = link)
    },
    error = function(cond) {
      return(cond)
    }
  )
  if (inherits(fit, "simpleError")) {
    stop("gee failed with error: ", fit[["message"]])
  }

  #---CHECKS---

  if (fit$cond == FALSE) {
    stop("standardize_gee is only implemented for gee object with cond=TRUE. For cond=FALSE, use stdGlm.")
  }
  link <- summary(fit)$link
  if (link != "identity" && link != "log") {
    stop("standardize_gee is only implemented for gee object with identity link or log link.")
  }

  #---PREPARATION---

  formula <- fit$formula
  weights <- rep(1, nrow(fit$x)) # gee does not allow for weights
  npar <- length(fit$coef)

  # Delete rows that did not contribute to the model fit,
  # e.g. missing data or not in subset for fit.
  m <- fit$x
  data <- data[match(rownames(m), rownames(data)), ]
  n <- nrow(data)

  ncluster <- length(unique(data[, clusterid]))

  nX <- nrow(valuesout)

  # Check if model.matrix works with object=formula. If formula contains splines,
  # then neither object=formula nor object=fit will not work when no variation
  # in exposure, since model.matrix needs to retrieve Boundary.knots
  # from terms(fit). Then fit glm so can use model.matrix with object=terms(fit.glm).
  data.x <- data
  data.x[, exposure_names] <- exposure[min(which(!is.na(exposure)))]
  m.x <- try(expr = model.matrix(object = formula, data = data.x), silent = TRUE)
  contains.splines <- FALSE
  if (!is.matrix(m.x)) {
    contains.splines <- TRUE
    environment(formula) <- new.env()
    fit.glm <- glm(formula = formula, data = data)
  }

  #---ESTIMATES OF INTERCEPTS---

  h <- as.vector(fit$x %*% matrix(fit$coef, nrow = npar, ncol = 1))
  dh.dbeta <- fit$x
  if (link == "identity") {
    r <- fit$y - h
    a <- ave(x = r, data[, clusterid], FUN = mean)
  }
  if (link == "log") {
    r <- fit$y * exp(-h)
    a <- log(ave(x = r, data[, clusterid], FUN = mean))
  }

  #---ESTIMATES OF MEANS AT VALUES SPECIFIED BY x ---

  pred <- matrix(nrow = n, ncol = nX)
  SI.beta <- matrix(nrow = nX, ncol = npar)
  for (i in 1:nX) {
    data.x <- do.call("transform", c(
      list(data),
      valuesout[i, , drop = FALSE]
    ))
    if (contains.splines) {
      m.x <- model.matrix(object = terms(fit.glm), data = data.x)[, -1, drop = FALSE]
    } else {
      m.x <- model.matrix(object = formula, data = data.x)[, -1, drop = FALSE]
    }
    h.x <- as.vector(m.x %*% matrix(fit$coef, nrow = npar, ncol = 1))
    dh.x.dbeta <- m.x
    eta <- a + h.x
    if (link == "identity") {
      mu <- eta
      dmu.deta <- rep(1, n)
      da.dbeta <- -apply(X = dh.dbeta, MARGIN = 2, FUN = ave, data[, clusterid])
    }
    if (link == "log") {
      mu <- exp(eta)
      dmu.deta <- mu
      da.dbeta <- -apply(X = r * dh.dbeta, MARGIN = 2, FUN = ave, data[, clusterid]) /
        exp(a)
    }
    pred[, i] <- mu
    deta.dbeta <- da.dbeta + dh.x.dbeta
    # When link=="log", exp(a) will be 0 if y=0 for all subjects in the cluster.
    # This causes da.dbeta and deta.dbeta to be NA, but they should be 0.
    deta.dbeta[is.na(deta.dbeta)] <- 0
    dmu.dbeta <- dmu.deta * deta.dbeta
    SI.beta[i, ] <- colMeans(weights * dmu.dbeta)
  }
  est <- colSums(weights * pred, na.rm = TRUE) /
    sum(weights)

  #---VARIANCE OF MEANS AT VALUES SPECIFIED BY x---

  ores <- weights * fit$x * fit$res
  mres <- weights * (pred - matrix(rep(est, each = n), nrow = n, ncol = nX))
  res <- cbind(mres, ores)
  res <- aggr(x = res, clusters = data[, clusterid])
  J <- var(res, na.rm = TRUE)

  SI <- cbind(-diag(nX) * mean(weights), SI.beta)
  oI <- cbind(
    matrix(0, nrow = npar, ncol = nX),
    -t(fit$x) %*% (weights * fit$d.res) / n
  )
  I <- rbind(SI, oI)
  V <- (solve(I) %*% J %*% t(solve(I)) * ncluster / n^2)[1:nX, 1:nX]
  vcov <- V

  fit_exposure <- contrast <- reference <- NULL
  variance <- vcov
  estimates <- est
  fit_outcome <- fit

  ## Add names to asymptotic covariance matrix
  rownames(variance) <- colnames(variance) <-
    do.call("paste", c(lapply(seq_len(ncol(valuesout)), function(i) {
      paste(names(valuesout)[i], "=", round(valuesout[[i]], 2L))
    }), sep = "; "))

  valuesout[["se"]] <- sqrt(diag(variance))
  valuesout[["estimates"]] <- estimates
  res <- list(
    estimates = valuesout,
    covariance = variance,
    fit_outcome = fit_outcome,
    fit_exposure = fit_exposure,
    exposure_names = exposure_names
  )

  format_result_standardize(
    res,
    contrasts,
    reference,
    transforms,
    ci_type,
    ci_level,
    "std_glm",
    "summary_std_glm"
  )
}
