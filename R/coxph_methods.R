#' @title Regression standardization in Cox proportional hazards models
#'
#' @description \code{stdCoxph} performs regression standardization in Cox proportional
#' hazards models, at specified values of the exposure, over the sample
#' covariate distribution. Let \eqn{T}, \eqn{X}, and \eqn{Z} be the survival
#' outcome, the exposure, and a vector of covariates, respectively.
#' \code{stdCoxph} uses a fitted Cox proportional hazards model to estimate the
#' standardized survival function \eqn{\theta(t,x)=E\{S(t|X=x,Z)\}}, where
#' \eqn{t} is a specific value of \eqn{T}, \eqn{x} is a specific value of
#' \eqn{X}, and the expectation is over the marginal distribution of \eqn{Z}.
#'
#' @details \code{stdCoxph} assumes that a Cox proportional hazards model
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
#'
#' @param fit an object of class \code{"coxph"}, as returned by the
#' \code{coxph} function in the \pkg{survival} package, but without special
#' terms \code{strata}, \code{cluster} or \code{tt}.  Only \code{breslow}
#' method for handling ties is allowed. If arguments \code{weights} and/or
#' \code{subset} are used when fitting the model, then the same weights and
#' subset are used in \code{stdGlm}.
#' @param data a data frame containing the variables in the model. This should
#' be the same data frame as was used to fit the model in \code{fit}.
#' @param X a string containing the name of the exposure variable \eqn{X} in
#' \code{data}.
#' @param x an optional vector containing the specific values of \eqn{X} at
#' which to estimate the standardized survival function. If \eqn{X} is binary
#' (0/1) or a factor, then \code{x} defaults to all values of \eqn{X}. If
#' \eqn{X} is numeric, then \code{x} defaults to the mean of \eqn{X}. If
#' \code{x} is set to \code{NA}, then \eqn{X} is not altered. This produces an
#' estimate of the marginal survival function \eqn{S(t)=E\{S(t|X,Z)\}}.
#' @param t an optional vector containing the specific values of \eqn{T} at
#' which to estimate the standardized survival function. It defaults to all the
#' observed event times in \code{data}.
#' @param clusterid an optional string containing the name of a cluster
#' identification variable when data are clustered.
#' @param subsetnew an optional logical statement specifying a subset of
#' observations to be used in the standardization. This set is assumed to be a
#' subset of the subset (if any) that was used to fit the regression model.
#' @return An object of class \code{"stdCoxph"} is a list containing
#' \item{call}{ the matched call.  } \item{input}{ \code{input} is a list
#' containing all input arguments.  } \item{est}{ a matrix with
#' \code{length(t)} rows and \code{length(x)} columns, where the element on row
#' \code{i} and column \code{j} is equal to
#' \eqn{\hat{\theta}}(\code{t[i],x[j]}).  } \item{vcov}{ a list with
#' \code{length(t)} elements. Each element is a square matrix with
#' \code{length(x)} rows. In the \code{k:th} matrix, the element on row
#' \code{i} and column \code{j} is the (estimated) covariance of
#' \eqn{\hat{\theta}}(\code{t[k]},\code{x[i]}) and
#' \eqn{\hat{\theta}}(\code{t[k]},\code{x[j]}).  }
#' @note Standardized survival functions are sometimes referred to as (direct)
#' adjusted survival functions in the literature.
#'
#' \code{stdCoxph} does not currently handle time-varying exposures or
#' covariates.
#'
#' \code{stdCoxph} internally loops over all values in the \code{t} argument.
#' Therefore, the function will usually be considerably faster if
#' \code{length(t)} is small.
#'
#' The variance calculation performed by \code{stdCoxph} does not condition on
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
#' X <- rnorm(n, mean=Z)
#' T <- rexp(n, rate=exp(X+Z+X*Z)) #survival time
#' C <- rexp(n, rate=exp(X+Z+X*Z)) #censoring time
#' U <- pmin(T, C) #time at risk
#' D <- as.numeric(T < C) #event indicator
#' dd <- data.frame(Z, X, U, D)
#' fit <- coxph(formula=Surv(U, D)~X+Z+X*Z, data=dd, method="breslow")
#' fit.std <- stdCoxph(fit=fit, data=dd, X="X", x=seq(-1,1,0.5), t=1:5)
#' print(summary(fit.std, t=3))
#' plot(fit.std)
#'
#' @export stdCoxph
stdCoxph <- function(fit, data, X, x, t, clusterid, subsetnew){
  call <- match.call()

  #---PREPARATION---
  if(!fit$method=="breslow")
    stop("Only breslow method for handling ties is allowed.", call.=FALSE)
  specials <- pmatch(c("strata(","cluster(","tt("), attr(terms(fit$formula),
                                                         "variables"))
  if(any(!is.na(specials)))
    stop("No special terms are allowed in the formula")

  formula <- fit$formula
  npar <- length(fit$coef)
  fit.detail <- coxph.detail(object=fit)

  #Delete rows that did not contribute to the model fit,
  #e.g. missing data or not in subset for fit.
  #Need to have object=fit in model.matrix, since neither object=formula nor
  #object=terms(fit) will remove rows not in subset.
  m <- model.matrix(object=fit)
  data <- data[match(rownames(m), rownames(data)), ]
  n <- nrow(data)

  #Make new subset if supplied.
  subsetnew <-
    if(missing(subsetnew))
      rep(1, n)
  else
    as.numeric(eval(substitute(subsetnew), data, parent.frame()))
  input <- as.list(environment())

  if(is.null(fit$weights))
    weights <- rep(1, nrow(data))
  else
    weights <- fit$weights

  #Can write code more generally with
  #if(missing(clusters)) clusters <- 1:nrow(data)
  #but a problem when constructing meat in sandwich formula:
  #must always aggregate, which is a bit slow, even though much faster
  #when using data.table than the aggregate function.
  if(!missing(clusterid))
    ncluster <- length(unique(data[, clusterid]))

  #Assign values to x and reference if not supplied.
  #Make sure x is a factor if data[, X] is a factor
  #with the same levels as data[, X].
  if(missing(x)){
    if(is.factor(data[, X]))
      x <- as.factor(levels(data[, X]))
    if(is.numeric(data[, X]))
      if(is.binary(data[, X]))
        x <- c(0, 1)
    else
      x <- round(mean(data[, X], na.rm=TRUE), 2)
  }
  else{
    if(is.factor(x)){
      temp <- x
      levels(x) <- levels(data[, X])
      x[1:length(x)] <- temp
    }
    else{
      if(is.factor(data[, X])){
        x <- factor(x)
        temp <- x
        levels(x) <- levels(data[, X])
        x[1:length(x)] <- temp
      }
    }
  }
  input$x <- x
  nX <- length(x)

  #Assign value to t if missing.
  if(missing(t))
    t <- fit.detail$time
  input$t <- t
  nt <- length(t)

  if(sum(fit.detail$time<=min(t))==0)
    stop("No events before first value in t", call.=FALSE)

  est <- matrix(nrow=nt, ncol=nX)
  vcov <- vector(mode="list", length=nt)
  H <- Hfun(fit=fit, data=data, fit.detail=fit.detail)
  sandwich.fit <- sandwich(fit=fit, data=data, weights=weights, t=t,
                           fit.detail=fit.detail)

  #---LOOP OVER nt

  for(j in 1:nt){

    if(t[j]==0){
      est[j, ] <- 1
      vcov[[j]] <- matrix(0, nrow=nX, ncol=nX)
    }
    else{

      #---ESTIMATES OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---

      si <- matrix(nrow=n, ncol=nX)
      PredX <- matrix(nrow=n, ncol=nX)
      tempmat <- matrix(nrow=nX, ncol=npar)
      for(i in 1:nX){
        data.x <- data
        if(!is.na(x[i]))
          data.x[, X] <- x[i]
        predX <- predict(object=fit, newdata=data.x, type="risk")
        si[, i] <- exp(-H(t[j])*predX)
        PredX[, i] <- predX
        #Need terms(fit) here. If formula contains splines,
        #then fit or formula will not work when no variation in the exposure,
        #since model.matrix need to retrieve Boundary.knots from terms(fit).
        #Also need to center, since everything else is centered.
        #Note: need to center at factual X, not counterfactual x,
        #since baseline is computed at mean of factual X.
        m.x <- model.matrix(object=terms(fit), data=data.x)[, -1, drop=FALSE]
        m <- model.matrix(object=terms(fit), data=data)[, -1, drop=FALSE]
        m <- matrix(colMeans(m), nrow=nrow(m), ncol=ncol(m), byrow=TRUE)
        m.x <- m.x-m
        tempmat[i, ] <- colMeans(m.x*predX*si[, i]*subsetnew*weights)
      }
      est[j, ] <- colSums(subsetnew*weights*si, na.rm=TRUE)/
        sum(subsetnew*weights)

      #---VARIANCE OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x, ---

      sres <- subsetnew*weights*(si-matrix(rep(est[j, ], each=n), nrow=n, ncol=nX))
      ores <- sandwich.fit$U[, c(1:npar, npar+j)]
      res <- cbind(sres, ores)
      if(!missing(clusterid))
        res <- aggr(x=res, clusters=data[, clusterid])
      J <- var(res, na.rm=TRUE)
      SI <- cbind(-diag(nX)*mean(subsetnew*weights), -tempmat*H(t[j]),
                  -colMeans(PredX*si*subsetnew*weights))
      #This is why the user cannot use term cluster; then -solve(vcov(object=fit))/n
      #will not be the bread in the sandwich.
      oI <- cbind(matrix(0, nrow=npar+1, ncol=nX),
                  sandwich.fit$I[c(1:npar, npar+j), c(1:npar, npar+j)])
      I <- rbind(SI, oI)

      if(missing(clusterid))
        V <- (solve(I)%*%J%*%t(solve(I))/n)[1:nX, 1:nX]
      else
        V <- (solve(I)%*%J%*%t(solve(I))*ncluster/n^2)[1:nX, 1:nX]
      vcov[[j]] <- V

    }

  }

  out <- list(call=call, input=input, est=est, vcov=vcov)

  #---OUTPUT---

  class(out) <- "stdCoxph"
  return(out)
}

#' @title Summarizes Cox regression standardization fit
#'
#' @description This is a \code{summary} method for class \code{"stdCoxph"}.
#'
#' @param object an object of class \code{"stdCoxph"}.
#' @param t numeric, indicating the times at which to summarize. It defaults to
#' the specified value(s) of the argument \code{t} in the \code{stdCox}
#' function.
#' @param CI.type string, indicating the type of confidence intervals. Either
#' "plain", which gives untransformed intervals, or "log", which gives
#' log-transformed intervals.
#' @param CI.level desired coverage probability of confidence intervals, on
#' decimal form.
#' @param transform a string. If set to \code{"log"}, \code{"logit"}, or
#' \code{"odds"}, the standardized survival function \eqn{\theta(t,x)} is
#' transformed into \eqn{\psi(t,x)=log\{\theta(t,x)\}},
#' \eqn{\psi(t,x)=log[\theta(t,x)/\{1-\theta(t,x)\}]}, or
#' \eqn{\psi(t,x)=\theta(t,x)/\{1-\theta(t,x)\}}, respectively. If left
#' unspecified, \eqn{\psi(t,x)=\theta(t,x)}.
#' @param contrast a string. If set to \code{"difference"} or \code{"ratio"},
#' then \eqn{\psi(t,x)-\psi(t,x_0)} or \eqn{\psi(t,x) / \psi(t,x_0)} are
#' constructed, where \eqn{x_0} is a reference level specified by the
#' \code{reference} argument.
#' @param reference must be specified if \code{contrast} is specified.
#' @param \dots not used.
#' @author Arvid Sjolander
#' @seealso \code{\link{stdCoxph}}
#' @examples
#'
#' ##See documentation for stdCoxph
#'
#' @rdname summary
#' @export summary.stdCoxph
#' @export
summary.stdCoxph <- function(object, t, CI.type="plain", CI.level=0.95,
                             transform=NULL, contrast=NULL, reference=NULL, ...){

  est.all <- object$est
  V.all <- object$vcov
  nX <- length(object$input$x)
  if(missing(t))
    t <- object$input$t
  nt <- length(t)

  est.table <- vector(mode="list", length=nt)
  for(j in 1:nt){

    if(min(abs(t[j]-object$input$t)) > sqrt(.Machine$double.eps))
      stop("The standardized survival function is not estimated at t",
           call.=FALSE)
    else
      k <- which.min(abs(t[j]-object$input$t))

    est <- est.all[k, ]
    V <- as.matrix(V.all[[k]])

    if(!is.null(transform)){
      if(transform=="log"){
        dtransform.dm <- diag(1/est, nrow=nX, ncol=nX)
        est <- log(est)
      }
      if(transform=="logit"){
        dtransform.dm <- diag(1/(est*(1-est)), nrow=nX, ncol=nX)
        est <- logit(est)
      }
      if(transform=="odds"){
        dtransform.dm <- diag(1/(1-est)^2, nrow=nX, ncol=nX)
        est <- odds(est)
      }
      V <- t(dtransform.dm)%*%V%*%dtransform.dm
    }

    if(!is.null(contrast)){
      if(is.null(reference))
        stop("When specifying contrast, reference must be specified as well")
      referencepos <- match(reference, object$input$x)
      if(is.na(referencepos))
        stop("reference must be a value in x")
      if(contrast=="difference"){
        dcontrast.dtransform <- diag(nX)
        dcontrast.dtransform[referencepos, ] <- -1
        dcontrast.dtransform[referencepos, referencepos] <- 0
        est <- est-est[referencepos]
      }
      if(contrast=="ratio"){
        dcontrast.dtransform <- diag(1/est[referencepos], nrow=nX, ncol=nX)
        dcontrast.dtransform[referencepos, ] <- -est/est[referencepos]^2
        dcontrast.dtransform[referencepos, referencepos] <- 1
        est <- est/est[referencepos]
      }
      V <- t(dcontrast.dtransform)%*%V%*%dcontrast.dtransform
      V[referencepos, ] <- 0
      V[, referencepos] <- 0
    }

    var <- diag(V)
    se <-  sqrt(var)
    conf.int <- CI(est=est, var=var, CI.type=CI.type, CI.level=CI.level)

    temp <- as.matrix(cbind(est, se, conf.int), nrow=length(est), ncol=4)
    dimnames(temp) <- list(object$input$x,
                           c("Estimate", "Std. Error", paste("lower",CI.level),
                             paste("upper", CI.level)))
    est.table[[j]] <- temp

  }
  if(is.factor(reference))
    reference <- as.character(reference)
  out <- c(object,
           list(est.table=est.table, tsum=t, transform=transform, contrast=contrast,
                reference=reference))
  class(out) <- "summary.stdCoxph"
  return(out)
}

#' @rdname print
#' @export print.stdCoxph
#' @export
print.stdCoxph <- function(x, ...){
  print(summary(x))
}

#' @title Prints summary of Cox regression standardization fit
#'
#' @description This is a \code{print} method for class \code{"summary.stdCoxph"}.
#'
#' @param x an object of class \code{"summary.stdCoxph"}.
#' @param \dots not used.
#' @author Arvid Sjolander
#' @seealso \code{\link{stdCoxph}}
#' @examples
#'
#'
#' ##See documentation for stdCoxph
#'
#' @rdname print
#' @export print.summary.stdCoxph
#' @export
print.summary.stdCoxph <- function(x, ...){
  nt <- length(x$tsum)
  for(j in 1:nt){
    cat("\nFormula: ")
    print(x$input$fit$formula, showEnv=FALSE)
    cat("Exposure: ", x$input$X, "\n")

    if(!is.null(x$transform))
      cat("Transform: ", x$transform,  "\n")
    if(!is.null(x$contrast)){
      cat("Reference level: ", x$input$X, "=", x$reference,  "\n")
      cat("Contrast: ", x$contrast,  "\n")
    }
    cat("Survival functions evaluated at t =", x$tsum[j], "\n")
    cat("\n")
    print(x$est.table[[j]], digits=3)
    cat("\n")
  }
}

#' @title Plots Cox regression standardization fit
#'
#' @description This is a \code{plot} method for class \code{"stdCoxph"}.
#'
#' @param x an object of class \code{"stdCoxph"}.
#' @param plot.CI logical, indicating whether confidence intervals should be
#' added to the plot.
#' @param CI.type string, indicating the type of confidence intervals. Either
#' "plain", which gives untransformed intervals, or "log", which gives
#' log-transformed intervals.
#' @param CI.level desired coverage probability of confidence intervals, on
#' decimal form.
#' @param transform a string. If set to \code{"log"}, \code{"logit"}, or
#' \code{"odds"}, the standardized survival function \eqn{\theta(t,x)} is
#' transformed into \eqn{\psi(t,x)=log\{\theta(t,x)\}},
#' \eqn{\psi(t,x)=log[\theta(t,x)/\{1-\theta(t,x)\}]}, or
#' \eqn{\psi(t,x)=\theta(t,x)/\{1-\theta(t,x)\}}, respectively. If left
#' unspecified, \eqn{\psi(t,x)=\theta(t,x)}.
#' @param contrast a string. If set to \code{"difference"} or \code{"ratio"},
#' then \eqn{\psi(t,x)-\psi(t,x_0)} or \eqn{\psi(t,x) / \psi(t,x_0)} are
#' constructed, where \eqn{x_0} is a reference level specified by the
#' \code{reference} argument.
#' @param reference must be specified if \code{contrast} is specified.
#' @param legendpos position of the legend; see help for \code{legend}.
#' @param \dots further arguments passed on to plot.default.
#' @author Arvid Sjolander
#' @seealso \code{\link{stdCoxph}}
#' @examples
#'
#' ##See documentation for stdCoxph
#'
#' @rdname plot
#' @export plot.stdCoxph
#' @export
plot.stdCoxph <- function(x, plot.CI=TRUE, CI.type="plain", CI.level=0.95,
                          transform=NULL, contrast=NULL, reference=NULL, legendpos="bottomleft", ...){

  object <- x
  x <- object$input$x

  dots <- list(...)

  xlab <- "t"

  if(is.factor(reference))
    reference <- as.character(reference)

  if(is.null(contrast)){
    if(is.null(transform))
      ylab <- expression(S(t))
    else{
      if(transform=="log")
        ylab <- expression(paste(log, "{", S(t), "}", sep=""))
      if(transform=="logit")
        ylab <- expression(paste(logit, "{", S(t), "}", sep=""))
      if(transform=="odds")
        ylab <- expression(paste(S(t), "/{", 1-S(t), "}", sep=""))
    }
  }
  else{
    if(contrast=="difference"){
      if(is.null(transform))
        ylab <- c(bquote(paste(S(t), "-", S[.(reference)](t))), expression())
      else{
        if(transform=="log")
          ylab <- c(bquote(paste(log, "{", S(t), "}-", log, "{",
                                 S[.(reference)](t), "}", sep="")), expression())
        if(transform=="logit")
          ylab <- c(bquote(paste(logit, "{", S(t), "}-", logit,
                                 "{", S[.(reference)](t), "}", sep="")), expression())
        if(transform=="odds")
          ylab <- c(bquote(paste(S(t), "/{", 1-S(t), "}-",
                                 S[.(reference)](t), "/{", 1-S[.(reference)](t),
                                 "}", sep="")), expression())
      }
    }
    if(contrast=="ratio"){
      if(is.null(transform))
        ylab <- c(bquote(paste(S(t), " / ", S[.(reference)](t), sep="")),
                  expression())
      else{
        if(transform=="log")
          ylab <- c(bquote(paste(log, "{", S(t), "} / ", log,
                                 "{", S[.(reference)](t), "}", sep="")), expression())
        if(transform=="logit")
          ylab <- c(bquote(paste(logit, "{", S(t), "} / ", logit,
                                 "{", S[.(reference)](t), "}", sep="")), expression())
        if(transform=="odds")
          ylab <- c(bquote(paste("[", S(t), "/{", 1-S(t), "}] / [",
                                 S[.(reference)](t), "/{", 1-S[.(reference)](t),
                                 "}]", sep="")), expression())
      }
    }
  }

  t <- object$input$t
  nt <- length(t)
  nX <- length(x)

  sum.obj <- summary(object=object, CI.type=CI.type, CI.level=CI.level,
                     transform=transform, contrast=contrast, reference=reference)
  temp <- Reduce(f=rbind, x=sum.obj$est.table)
  est <- matrix(temp[, 1], nrow=nt, ncol=nX, byrow=TRUE)
  lower <- matrix(temp[, 3], nrow=nt, ncol=nX, byrow=TRUE)
  upper <- matrix(temp[, 4], nrow=nt, ncol=nX, byrow=TRUE)

  if(plot.CI)
    ylim <- c(min(lower), max(upper))
  else
    ylim <- c(min(est), max(est))
  args <- list(x=object$input$t, y=rep(0, length(t)), xlab=xlab, ylab=ylab,
               ylim=ylim, type="n")
  args[names(dots)] <- dots
  do.call("plot", args=args)
  legend <- NULL
  for(i in 1:nX){
    lines(t, est[, i], col=i)
    if(plot.CI){
      lines(t, upper[, i], lty="dashed", col=i)
      lines(t, lower[, i], lty="dashed", col=i)
    }
    temp <- as.character(x[i])
    legend <- c(legend, paste(object$input$X, "=", object$input$x[i]))
  }
  if(is.na(match("ylim",names(args))))
    yl <- ylim[2]
  else
    yl <- args$ylim[2]
  legend(x=legendpos, legend=legend, lty=rep(1, length(x)), col=1:length(x),
         bty="n")

}
