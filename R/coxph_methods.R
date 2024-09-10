#' @title Regression standardization in Cox proportional hazards models
#'
#' @description \code{standardize_coxph} performs regression standardization in Cox proportional
#' hazards models at specified values of the exposure over the sample
#' covariate distribution. Let \eqn{T}, \eqn{X}, and \eqn{Z} be the survival
#' outcome, the exposure, and a vector of covariates, respectively.
#' \code{standardize_coxph} fits a Cox proportional hazards model and the Breslow estimator
#' of the baseline hazard in order to estimate the
#' standardized survival function \eqn{\theta(t,x)=E\{S(t|X=x,Z)\}} when \code{measure = "survival"} or the standardized restricted mean survival up to time t \eqn{\theta(t, x) = E\{\int_0^t S(u|X = x, Z) du\}} when \code{measure = "rmean"}, where
#' \eqn{t} is a specific value of \eqn{T}, \eqn{x} is a specific value of
#' \eqn{X}, and the expectation is over the marginal distribution of \eqn{Z}.
#'
#' @details \code{standardize_coxph} fits the Cox proportional hazards model
#' \deqn{\lambda(t|X,Z)=\lambda_0(t)\exp\{h(X,Z;\beta)\}.}
#' Breslow's estimator of the cumulative baseline hazard
#' \eqn{\Lambda_0(t)=\int_0^t\lambda_0(u)du} is used together with the partial
#' likelihood estimate of \eqn{\beta} to obtain estimates of the survival
#' function \eqn{S(t|X=x,Z)} if \code{measure = "survival"}:
#' \deqn{\hat{S}(t|X=x,Z)=\exp[-\hat{\Lambda}_0(t)\exp\{h(X=x,Z;\hat{\beta})\}].}
#' For each \eqn{t} in the \code{t} argument and for each \eqn{x} in the
#' \code{x} argument, these estimates are averaged across all subjects (i.e.
#' all observed values of \eqn{Z}) to produce estimates
#' \deqn{\hat{\theta}(t,x)=\sum_{i=1}^n \hat{S}(t|X=x,Z_i)/n,} where \eqn{Z_i}
#' is the value of \eqn{Z} for subject \eqn{i}, \eqn{i=1,...,n}.  The variance
#' for \eqn{\hat{\theta}(t,x)} is obtained by the sandwich formula.
#'
#' If \code{measure = "rmean"}, then \eqn{\Lambda_0(t)=\int_0^t\lambda_0(u)du}
#' is used together with the partial
#' likelihood estimate of \eqn{\beta} to obtain estimates of the restricted mean survival
#' up to time t: \eqn{\int_0^t S(u|X=x,Z) du} for each element of \code{times}. The estimation
#' and inference is done using the method described in Chen and Tsiatis 2001.
#' Currently, we can only estimate the difference in RMST for a single binary
#' exposure. Two separate Cox models are fit for each level of the exposure,
#' which is expected to be coded as 0/1.
#' @return An object of class \code{std_surv}.
#' This is basically a list with components estimates and covariance estimates in \code{res}
#' Results for transformations, contrasts, references are stored in \code{res_contrasts}.
#' The output contains estimates for contrasts and confidence intervals for all
#' combinations of transforms and reference levels.
#' Obtain numeric results in a data frame with the \link{tidy} function.
#' @inherit standardize_glm
#' @param times A vector containing the specific values of \eqn{T} at
#' which to estimate the standardized survival function.
#' @param measure Either "survival" to estimate the survival function at times
#' or "rmean" for the restricted mean survival up to the largest of times.
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
#' @author Arvid Sjölander, Adam Brand, Michael Sachs
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
#' Sjölander A. (2016). Regression standardization with the R-package stdReg.
#' \emph{European Journal of Epidemiology} \bold{31}(6), 563-574.
#'
#' Sjölander A. (2018). Estimation of causal effect measures with the R-package
#' stdReg. \emph{European Journal of Epidemiology} \bold{33}(9), 847-858.
#'
#' Chen, P. Y., Tsiatis, A. A. (2001). Causal inference on the difference of the restricted mean lifetime between two groups. \emph{Biometrics}, \bold{57}(4), 1030-1038.
#' @examples
#'
#'
#' require(survival)
#' set.seed(7)
#' n <- 300
#' Z <- rnorm(n)
#' Zbin <- rbinom(n, 1, .3)
#' X <- rnorm(n, mean = Z)
#' T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
#' C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
#' fact <- factor(sample(letters[1:3], n, replace = TRUE))
#' U <- pmin(T, C) # time at risk
#' D <- as.numeric(T < C) # event indicator
#' dd <- data.frame(Z, Zbin, X, U, D, fact)
#' fit.std.surv <- standardize_coxph(
#'   formula = Surv(U, D) ~ X + Z + X * Z,
#'   data = dd,
#'   values = list(X = seq(-1, 1, 0.5)),
#'   times = 1:5
#' )
#' print(fit.std.surv)
#' plot(fit.std.surv)
#' tidy(fit.std.surv)
#'
#' fit.std <- standardize_coxph(
#'   formula = Surv(U, D) ~ X + Zbin + X * Zbin + fact,
#'   data = dd,
#'   values = list(Zbin = 0:1),
#'   times = 1.5,
#'   measure = "rmean",
#'   contrast = "difference",
#'   reference = 0
#' )
#' print(fit.std)
#' tidy(fit.std)
#'
#' @export standardize_coxph
standardize_coxph <- function(formula,
                              data,
                              values,
                              times,
                              measure = c("survival", "rmean"),
                              clusterid,
                              ci_level = 0.95,
                              ci_type = "plain",
                              contrasts = NULL,
                              family = "gaussian",
                              reference = NULL,
                              transforms = NULL) {
  call <- match.call()

  measure <- match.arg(measure)
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

  if(measure == "rmean") {
    if (length(values) > 1L) {
      stop("there has to be only one exposure
              with the restricted mean survival estimator")
    }

    ## Check that exposure is binary
    if (!is.binary(data[[names(values)]]) || nrow(valuesout) != 2L) {
      stop("the exposure has to be binary (0 or 1)")
    }

    #---PREPARATION---
    specials <- pmatch(c("strata(", "cluster(", "tt("), attr(
      terms(formula),
      "variables"
    ))
    if (any(!is.na(specials))) {
      stop("No special terms are allowed in the formula")
    }
    tstar <- max(times)

    if(length(times) > 1L) {

      warning("Can only compute restricted mean for one time, using maximum supplied t=", max(times))
    }

    expname <- exposure_names <- names(values)[1]
    newformula <- update(formula, reformulate(attr(drop.terms(terms(formula),
               grep(expname, attr(terms(formula), "term.labels"))), "term.labels")))

    data0 <- data[data[[expname]] == 0,]
    data1 <- data[data[[expname]] == 1,]

    coxfit0 <- tryCatch(coxph(newformula, data = data0, model = TRUE),
                        error = function(e) e)
    if (inherits(coxfit0, "simpleError")) {
      stop("Cox fit in group 0 failed with error: ", coxfit0[["message"]])
    }

    coxfit1 <- tryCatch(coxph(newformula, data = data1, model = TRUE),
                        error = function(e) e)
    if (inherits(coxfit1, "simpleError")) {
      stop("Cox fit in group 1 failed with error: ", coxfit1[["message"]])
    }

    m0 <- model.matrix(object = coxfit0, data = data0)
    m1 <- model.matrix(object = coxfit1, data = data1)
    data <- data[match(sort(c(rownames(m0), rownames(m1))), rownames(data)), ]
    n <- nrow(data)

    input <- as.list(environment())

    input$expname <- NULL
    mf <- model.frame(newformula, data)
    survobj <- mf[,1]

    dmat <- model.matrix(mf, data = data)[, -1, drop = FALSE]

    ncovs <- ncol(dmat)

    input$times <- tstar

    etimes <- sort(survobj[survobj[, 2] == 1,1])
    etimes <- etimes[etimes <= tstar]
    Ai <- data[[expname]][match(etimes, survobj[,1])]

    ## get the names provided to Surv
    ##
    if(is.null(attr(terms(newformula), "variables")[[2]] |> names())) {
      timename <- deparse(attr(terms(newformula), "variables")[[2]][[2]])
      eventname <- deparse(attr(terms(newformula), "variables")[[2]][[3]])
    } else {
      timename <- deparse(attr(terms(newformula), "variables")[[2]][["time"]])
      eventname <- deparse(attr(terms(newformula), "variables")[[2]][["event"]])
    }

    basedat <- do.call(rbind, lapply(etimes, \(tt) {
      ret <- mf[, -1, drop = FALSE]
      ret$id <- 1:nrow(ret)
      ret[[timename]] <- tt
      ret[[eventname]] <- 1
      ret
    }))

    covs <- names(mf)[-1]

    breslow0 <- matrix(predict(coxfit0, newdata = basedat, type = "expected"),
                       nrow = nrow(data), ncol = length(etimes))
    breslow1 <- matrix(predict(coxfit1, newdata = basedat, type = "expected"),
                       nrow = nrow(data), ncol = length(etimes))

    risk0 <- predict(coxfit0, newdata = data, type = "risk", reference = "zero")
    risk1 <- predict(coxfit1, newdata = data, type = "risk", reference = "zero")

    S0hat <- colMeans(exp(-breslow0))
    S1hat <- colMeans(exp(-breslow1))

    rmst0 <- rsum(c(1,S0hat), c(0,etimes), tstar)
    rmst1 <- rsum(c(1,S1hat), c(0,etimes), tstar)

    diffrmst <- rmst1 - rmst0

    ## individual counting process

    Nit <- do.call(cbind, lapply(etimes, \(tt){

      ifelse(survobj[,1] <= tt & survobj[,2] == 1, 1, 0)

    }))

    Yit <- outer(survobj[,1], etimes, \(x,y) 1.0 * (x >= y))


    h0hattmp <- matrix(diff(c(0, etimes, tstar)),
                       nrow = nrow(data), ncol = length(etimes) + 1, byrow = TRUE) *
      cbind(1,exp(-breslow0)) * matrix(risk0, nrow = nrow(data),
                                       ncol = length(etimes)+1)
    h0hat <- colMeans(rowcumSums(h0hattmp[, ncol(h0hattmp):1])[,ncol(h0hattmp):1])[-1]

    h1hattmp <- matrix(diff(c(0, etimes, tstar)), nrow = nrow(data),
                       ncol = length(etimes) + 1, byrow = TRUE) *
      cbind(1,exp(-breslow1)) * matrix(risk1, nrow = nrow(data),
                                       ncol = length(etimes)+1)
    h1hat <- colMeans(rowcumSums(h1hattmp[, ncol(h1hattmp):1])[, ncol(h1hattmp):1])[-1]


    SS0.t.0 <- colMeans(matrix(risk0 * (data[[expname]] == 0), nrow = nrow(data),
                               ncol = length(etimes), byrow = FALSE) * Yit)


    SS1.t.0 <- matrix(NA, nrow = length(etimes), ncol = ncovs)
    SS2.t.0 <- array(NA, dim = c(length(etimes), ncovs, ncovs))

    SS0.t.1 <- colMeans(matrix(risk1 * (data[[expname]] == 1), nrow = nrow(data),
                               ncol = length(etimes), byrow = FALSE) * Yit)


    SS1.t.1 <- matrix(NA, nrow = length(etimes), ncol = ncovs)
    SS2.t.1 <- array(NA, dim = c(length(etimes), ncovs, ncovs))

    Sigma0 <- Sigma1 <- matrix(0, nrow = ncovs, ncol = ncovs)


    ZZmat <-  do.call(rbind,
                      apply(dmat, MARGIN = 1, FUN = \(x) c(x %*% t(x)),
                            simplify = FALSE))


    for(i in 1:length(etimes)) {

      SS1.t.0[i, ] <- colMeans(matrix((data[[expname]] == 0) * risk0 * Yit[, i],
                                      nrow = nrow(data),
                                      ncol = ncovs) * dmat)



      SS2.t.0[i,,] <- matrix(colMeans(matrix((data[[expname]] == 0) * risk0 * Yit[, i],
                                             nrow = nrow(data),
                                             ncol = ncovs^2) *
                                        ZZmat),
                             nrow = ncovs,
                             ncol = ncovs)

      SS1.t.1[i, ] <- colMeans(matrix((data[[expname]] == 1) * risk1 * Yit[, i],
                                      nrow = nrow(data),
                                      ncol = ncovs) * dmat)


      SS2.t.1[i,,] <- matrix(colMeans(matrix((data[[expname]] == 1) * risk1 * Yit[, i],
                                             nrow = nrow(data),
                                             ncol = ncovs^2) *
                                        ZZmat),
                             nrow = ncovs,
                             ncol = ncovs)


      Sigma0 <- Sigma0 + (1 - Ai[i]) * (SS2.t.0[i,,] / SS0.t.0[i] -
                                          (SS1.t.0[i,] / SS0.t.0[i]) %*% t(SS1.t.0[i,] / SS0.t.0[i]))
      Sigma1 <- Sigma1 + (Ai[i]) * (SS2.t.1[i,,] / SS0.t.1[i] -
                                      (SS1.t.1[i,] / SS0.t.1[i]) %*% t(SS1.t.1[i,] / SS0.t.1[i]))

    }

    Sigma0 <- Sigma0 / sum(1 - data[[expname]])
    Sigma1 <- Sigma1 / sum(data[[expname]])


    Zbar <- SS1.t.0 / matrix(SS0.t.0, nrow = length(SS0.t.0), ncol = ncol(SS1.t.0))
    Zbar1 <- SS1.t.1 / matrix(SS0.t.1, nrow = length(SS0.t.1), ncol = ncol(SS1.t.1))

    denom.haz0 <- colSums(matrix(risk0 * (data[[expname]] == 0),
                                 nrow = length(risk0), ncol = ncol(Yit)) * Yit)
    denom.haz1 <- colSums(matrix(risk1 * (data[[expname]] == 1),
                                 nrow = length(risk1), ncol = ncol(Yit)) * Yit)
    g0i <- g1i <- matrix(NA, nrow = nrow(data), ncol = ncovs)


    for(i in 1:nrow(data)) {


      inmat <- (1 - Ai) * (Zbar - dmat[replicate(nrow(Zbar), i), ]) /
        matrix(denom.haz0, nrow = length(denom.haz0), ncol = ncol(Zbar))

      inmat1 <- Ai * (Zbar1 - dmat[replicate(nrow(Zbar1), i), ]) /
        matrix(denom.haz1, nrow = length(denom.haz1), ncol = ncol(Zbar1))

      Lcurv0 <-  rbind(0, matrix(exp(-breslow0[i, ]) * risk0[i], nrow = nrow(inmat),
                                 ncol = ncol(inmat)) *apply(inmat, MARGIN = 2, cumsum))

      Lcurv1 <-  rbind(0, matrix(exp(-breslow1[i, ]) * risk1[i], nrow = nrow(inmat1),
                                 ncol = ncol(inmat1)) *apply(inmat1, MARGIN = 2, cumsum))

      g0i[i,] <- apply(Lcurv0, MARGIN = 2, rsum, c(0, etimes), tstar)
      g1i[i,] <- apply(Lcurv1, MARGIN = 2, rsum, c(0, etimes), tstar)


    }

    g0 <- (solve(Sigma0) %*% colSums(g0i)) / sum(1 - data[[expname]])
    g1 <- (solve(Sigma1) %*% colSums(g1i)) / sum(data[[expname]])

    var0 <- (c((sum(1 - data[[expname]]) / nrow(data)) * t(g0) %*% Sigma0 %*% g0) +
      sum((1 - Ai) * (h0hat^2 / SS0.t.0) / denom.haz0) +
      (nrow(data) - 1) / (nrow(data)) * var(apply(cbind(1,exp(-breslow0)),
                                                  MARGIN = 1, rsum, c(0,etimes),tstar))) /
      nrow(data)

    var1 <- (c((sum(data[[expname]]) / nrow(data)) * t(g1) %*% Sigma1 %*% g1) +
                sum(Ai * (h1hat^2 / SS0.t.1) / denom.haz1) +
                (nrow(data) - 1) / (nrow(data)) * var(apply(cbind(1,exp(-breslow1)),
                                                            MARGIN = 1, rsum, c(0,etimes),tstar))) /
  nrow(data)

    covar <- (nrow(data) - 1) / (nrow(data)^2) * cov(apply(cbind(1,exp(-breslow1)),
                                                         MARGIN = 1, rsum, c(0,etimes),tstar),
                                                     apply(cbind(1,exp(-breslow0)),
                                                           MARGIN = 1, rsum, c(0,etimes),tstar))

    #var.diff <- var1 + var0 - 2 * covar

    est <- rbind(c(rmst0, rmst1))
    vcov <- list(matrix(c(var0, covar, covar, var1), nrow = 2, ncol = 2))

    }  else if(measure == "survival") {
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
  exposure_names <- names(values)[1]

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
    measure_name <- if(measure == "survival"){
      "survival function"
    } else if (measure == "rmean") {
        "restricted mean survival"
      }
    stop(sprintf("You have to specify the times (t) at which to estimate the standardized %s", measure_name))
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

  } else {
    stop("Measure ", measure, " not supported. Did you mean survival or rmean?")
  }

  out <- list(call = call, input = input, measure = measure,
              est = est, vcov = vcov)
  class(out) <- "std_coxph"

  #---OUTPUT---
  format_result_standardize(
    out,
    contrasts,
    reference,
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
  measure_name <- if(object$measure == "survival"){
    "survival function"
  } else if (object$measure == "rmean") {
    "restricted mean survival"
  }
  est.table <- vector(mode = "list", length = nt)
  for (j in 1:nt) {
    if (min(abs(times[j] - object$input$times)) > sqrt(.Machine$double.eps)) {
      stop(sprintf("The standardized %s is not estimated at times", measure_name),
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
        if(out$measure == "rmean"){
          stop("Transform logit not available for restricted mean survival.")
        }
        dtransform.dm <- diag(1 / (est * (1 - est)), nrow = nX, ncol = nX)
        est <- logit(est)
      }
      if (transform == "odds") {

        if(out$measure == "rmean"){
          stop("Transform odds not available for restricted mean survival.")
        }
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

    temp <- data.frame(object$input$valuesout, est, se, conf.int)
    colnames(temp) <-
      c(
        colnames(object$input$valuesout)[1], "Estimate", "Std.Error", paste0("lower.", ci_level),
        paste("upper", ci_level)
      )

    est.table[[j]] <- temp
  }
  if (is.factor(reference)) {
    reference <- as.character(reference)
  }
  out <- c(
    object,
    list(
      est_table = est.table, times = times, measure = object$measure,
      transform = transform, contrast = contrast,
      reference = reference, ci_type = ci_type, ci_level = ci_level
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
  invisible(x)
}

print_summary_std_coxph <- function(x, ...) {
  nt <- length(x$times)
  for (j in 1:nt) {
    cat("\nFormula: ")
    if(x$measure == "rmean") {
      print(x$input$newformula, showEnv = FALSE)
      cat(" Fit separately by exposure", "\n")
    } else {
      print(x$input$fit$formula, showEnv = FALSE)
    }

    cat("Exposure: ", x$input$exposure_names, "\n")

    if (!is.null(x$transform)) {
      cat("Transform: ", x$transform, "\n")
    }
    if (!is.null(x$contrast)) {
      cat("Reference level: ", x$input$X, "=", x$reference, "\n")
      cat("Contrast: ", x$contrast, "\n")
    }
    measure_name <- switch(x$measure, survival = "Survival functions",
                           rmean = "Restricted mean survival")
    cat(sprintf("%s evaluated at t =", measure_name), x$times[j], "\n")
    cat("\n")
    print(x$est_table[[j]], digits = 3)
    cat("\n")
  }
}

#' @title Plots regression standardization fit
#' @description This is a \code{plot} method for class \code{"std_surv"}.
#' @param x An object of class \code{"std_surv"}.
#' @param plot_ci if \code{TRUE}, add the confidence intervals to the plot.
#' @param ci_type A string, indicating the type of confidence intervals. Either "plain", which
#' gives untransformed intervals, or "log", which gives log-transformed intervals.
#' @param ci_level Coverage probability of confidence intervals.
#' @param transform  If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized
#' mean \eqn{\theta(x)} is transformed into \eqn{\psi(x)=\log\{\theta(x)\}},
#' \eqn{\psi(x)=\log[\theta(x)/\{1-\theta(x)\}]}, or
#' \eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively. If left unspecified,
#' \eqn{\psi(x)=\theta(x)}.
#' @param contrast If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(x)-\psi(x_0)}
#' or \eqn{\psi(x) / \psi(x_0)} are constructed, where \eqn{x_0} is a reference
#' level specified by the \code{reference} argument.
#' If not \code{NULL}, a doubly robust estimator of the standardized estimator is used.
#' @param reference If \code{contrast} is specified, the desired reference level.
#' @param summary_fun For internal use only. Do not change.
#' @param legendpos position of the legend; see \link[graphics]{legend}.
#' @param \dots Unused.
#' @rdname plot.std_surv
#' @export plot.std_surv
#' @export
#' @returns None. Creates a plot as a side effect
plot.std_surv <- function(x, plot_ci = TRUE, ci_type = "plain", ci_level = 0.95,
                          transform = NULL, contrast = NULL,
                          reference = NULL, legendpos = "bottomleft",
                          summary_fun = "summary_std_coxph", ...) {
  object <- x$res
  if (ncol(object$input$valuesout) != 1) {
    stop("multiple exposures")
  } else {
    x <- object$input$valuesout[, 1]
  }

  if(object$measure == "rmean") {
    stop("Cannot plot restricted mean")
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

  temp <- Reduce(f = rbind, x = sum.obj$est_table)
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



#' Provide tidy output from a std_surv object for use in downstream computations
#'
#' Tidy summarizes information about the components of the standardized model fit.
#' @param x An object of class std_surv
#' @param ... Not currently used
#'
#' @returns A data.frame
#' @examples
#' require(survival)
#' set.seed(8)
#' n <- 300
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' time <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
#' C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
#' U <- pmin(time, C) # time at risk
#' D <- as.numeric(time < C) # event indicator
#' dd <- data.frame(Z, X, U, D)
#' x <- standardize_coxph(
#'   formula = Surv(U, D) ~ X + Z + X * Z,
#'   data = dd, values = list(X = seq(-1, 1, 0.5)), times = c(2,3,4)
#' )
#'
#' tidy(x)
#'
#' @export
tidy.std_surv <- function(x, ...) {

  stopifnot(inherits(x, "std_surv"))

  res_list <- lapply(x$res_contrast, \(xl) {

    for(i in seq_along(xl$est_table)){
     xl$est_table[[i]]$time <- xl$input$times[i]
    }
    tmpres <- do.call(rbind, xl$est_table)
    colnames(tmpres) <- make.names(colnames(tmpres))
    tmpres$contrast <- if(is.null(xl$contrast)) "none" else xl$contrast
    tmpres$transform <- if(is.null(xl$transform)) "identity" else xl$transform
    tmpres$measure <- xl$measure
    tmpres

  })

  do.call("rbind", res_list)

}
