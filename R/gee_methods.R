#' @title Regression standardization in conditional generalized estimating equations
#'
#' @description \code{stdGee} performs regression standardization in linear and log-linear
#' fixed effects models, at specified values of the exposure, over the sample
#' covariate distribution. Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome,
#' the exposure, and a vector of covariates, respectively. It is assumed that
#' data are clustered with a cluster indicator \eqn{i}. \code{stdGee} uses
#' fitted fixed effects model, with cluster-specific intercept \eqn{a_i} (see
#' \code{details}), to estimate the standardized mean
#' \eqn{\theta(x)=E\{E(Y|i,X=x,Z)\}}, where \eqn{x} is a specific value of
#' \eqn{X}, and the outer expectation is over the marginal distribution of
#' \eqn{(a_i,Z)}.
#'
#' @details \code{stdGee} assumes that a fixed effects model
#' \deqn{\eta\{E(Y|i,X,Z)\}=a_i+h(X,Z;\beta)} has been fitted. The link
#' function \eqn{\eta} is assumed to be the identity link or the log link. The
#' conditional generalized estimating equation (CGGE) estimate of \eqn{\beta}
#' is used to obtain estimates of the cluster-specific means:
#' \deqn{\hat{a}_i=\sum_{j=1}^{n_i}r_{ij}/n_i,} where
#' \deqn{r_{ij}=Y_{ij}-h(X_{ij},Z_{ij};\hat{\beta})} if \eqn{\eta} is the
#' identity link, and \deqn{r_{ij}=Y_{ij}exp\{-h(X_{ij},Z_{ij};\hat{\beta})\}}
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
#' @param fit an object of class \code{"gee"}, with argument \code{cond =
#' TRUE}, as returned by the \code{gee} function in the \pkg{drgee} package. If
#' arguments \code{weights} and/or \code{subset} are used when fitting the
#' model, then the same weights and subset are used in \code{stdGee}.
#' @param data a data frame containing the variables in the model. This should
#' be the same data frame as was used to fit the model in \code{fit}.
#' @param X a string containing the name of the exposure variable \eqn{X} in
#' \code{data}.
#' @param x an optional vector containing the specific values of \eqn{X} at
#' which to estimate the standardized mean. If \eqn{X} is binary (0/1) or a
#' factor, then \code{x} defaults to all values of \eqn{X}. If \eqn{X} is
#' numeric, then \code{x} defaults to the mean of \eqn{X}. If \code{x} is set
#' to \code{NA}, then \eqn{X} is not altered. This produces an estimate of the
#' marginal mean \eqn{E(Y)=E\{E(Y|X,Z)\}}.
#' @param clusterid an mandatory string containing the name of a cluster
#' identification variable. Must be identical to the clusterid variable used in
#' the gee call.
#' @param subsetnew an optional logical statement specifying a subset of
#' observations to be used in the standardization. This set is assumed to be a
#' subset of the subset (if any) that was used to fit the regression model.
#' @return An object of class \code{"stdGee"} is a list containing \item{call}{
#' the matched call.  } \item{input}{ \code{input} is a list containing all
#' input arguments.  } \item{est}{ a vector with length equal to
#' \code{length(x)}, where element \code{j} is equal to
#' \eqn{\hat{\theta}}(\code{x[j]}).  } \item{vcov}{ a square matrix with
#' \code{length(x)} rows, where the element on row \code{i} and column \code{j}
#' is the (estimated) covariance of \eqn{\hat{\theta}}(\code{x[i]}) and
#' \eqn{\hat{\theta}}(\code{x[j]}).  }
#' @note The variance calculation performed by \code{stdGee} does not condition
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
#' @author Arvid Sjolander.
#' @references Goetgeluk S. and Vansteelandt S. (2008). Conditional generalized
#' estimating equations for the analysis of clustered and longitudinal data.
#' \emph{Biometrics} \bold{64}(3), 772-780.
#'
#' Martin R.S. (2017). Estimation of average marginal effects in multiplicative
#' unobserved effects panel models. \emph{Economics Letters} \bold{160}, 16-19.
#'
#' Sjolander A. (2019). Estimation of marginal causal effects in the presence
#' of confounding by cluster. \emph{Biostatistics} doi:
#' 10.1093/biostatistics/kxz054
#' @examples
#'
#' require(drgee)
#'
#' set.seed(4)
#' n <- 300
#' ni <- 2
#' id <- rep(1:n, each=ni)
#' ai <- rep(rnorm(n), each=ni)
#' Z <- rnorm(n*ni)
#' X <- rnorm(n*ni, mean=ai+Z)
#' Y <- rnorm(n*ni, mean=ai+X+Z+0.1*X^2)
#' dd <- data.frame(id, Z, X, Y)
#' fit <- gee(formula=Y~X+Z+I(X^2), data=dd, clusterid="id", link="identity",
#'   cond=TRUE)
#' fit.std <- stdGee(fit=fit, data=dd, X="X", x=seq(-3,3,0.5), clusterid="id")
#' print(summary(fit.std, contrast="difference", reference=2))
#' plot(fit.std)
#'
#' @export stdGee
stdGee <- function(fit, data, X, x, clusterid, subsetnew){
  call <- match.call()

  #---CHECKS---

  if(fit$cond==FALSE)
    stop("stdGee is only implemented for gee object with cond=TRUE. For cond=FALSE, use stdGlm.")
  link <- summary(fit)$link
  if(link!="identity" & link!="log")
    stop("stdGee is only implemented for gee object with identity link or log link.")

  #---PREPARATION---

  formula <- fit$formula
  weights <- rep(1, nrow(fit$x)) #gee does not allow for weights
  npar <- length(fit$coef)

  #Delete rows that did not contribute to the model fit,
  #e.g. missing data or not in subset for fit.
  m <- fit$x
  data <- data[match(rownames(m), rownames(data)), ]
  n <- nrow(data)

  #Make new subset if supplied.
  subsetnew <-
    if(missing(subsetnew))
      rep(1, n)
  else
    as.numeric(eval(substitute(subsetnew), data, parent.frame()))
  input <- as.list(environment())

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

  #Check if model.matrix works with object=formula. If formula contains splines,
  #then neither object=formula nor object=fit will not work when no variation
  #in exposure, since model.matrix needs to retrieve Boundary.knots
  #from terms(fit). Then fit glm so can use model.matrix with object=terms(fit.glm).
  data.x <- data
  data.x[, X] <- x[min(which(!is.na(x)))]
  m.x <- try(expr=model.matrix(object=formula, data=data.x), silent=TRUE)
  contains.splines <- FALSE
  if(!is.matrix(m.x)){
    contains.splines <- TRUE
    environment(formula) <- new.env()
    fit.glm <- glm(formula=formula, data=data)
  }

  #---ESTIMATES OF INTERCEPTS---

  h <- as.vector(fit$x%*%matrix(fit$coef, nrow=npar, ncol=1))
  dh.dbeta <- fit$x
  if(link=="identity"){
    r <- fit$y-h
    a <- ave(x=r, data[, clusterid], FUN=mean)
  }
  if(link=="log"){
    r <- fit$y*exp(-h)
    a <- log(ave(x=r, data[, clusterid], FUN=mean))
  }

  #---ESTIMATES OF MEANS AT VALUES SPECIFIED BY x ---

  pred <- matrix(nrow=n, ncol=nX)
  SI.beta <- matrix(nrow=nX, ncol=npar)
  for(i in 1:nX){
    data.x <- data
    if(!is.na(x[i]))
      data.x[, X] <- x[i]
    if(contains.splines){
      m.x <- model.matrix(object=terms(fit.glm), data=data.x)[, -1, drop=FALSE]
    }
    else{
      m.x <- model.matrix(object=formula, data=data.x)[, -1, drop=FALSE]
    }
    h.x <- as.vector(m.x%*%matrix(fit$coef, nrow=npar, ncol=1))
    dh.x.dbeta <- m.x
    eta <- a+h.x
    if(link=="identity"){
      mu <- eta
      dmu.deta <- rep(1, n)
      da.dbeta <- -apply(X=dh.dbeta, MARGIN=2, FUN=ave, data[, clusterid])
    }
    if(link=="log"){
      mu <- exp(eta)
      dmu.deta <- mu
      da.dbeta <- -apply(X=r*dh.dbeta, MARGIN=2, FUN=ave, data[, clusterid])/
        exp(a)
    }
    pred[, i] <- mu
    deta.dbeta <- da.dbeta+dh.x.dbeta
    #When link=="log", exp(a) will be 0 if y=0 for all subjects in the cluster.
    #This causes da.dbeta and deta.dbeta to be NA, but they should be 0.
    deta.dbeta[is.na(deta.dbeta)] <- 0
    dmu.dbeta <- dmu.deta*deta.dbeta
    SI.beta[i, ] <- colMeans(subsetnew*weights*dmu.dbeta)
  }
  est <- colSums(subsetnew*weights*pred, na.rm=TRUE)/
    sum(subsetnew*weights)

  #---VARIANCE OF MEANS AT VALUES SPECIFIED BY x---

  ores <- weights*fit$x*fit$res
  mres <- subsetnew*weights*(pred-matrix(rep(est, each=n), nrow=n, ncol=nX))
  res <- cbind(mres, ores)
  res <- aggr(x=res, clusters=data[, clusterid])
  J <- var(res, na.rm=TRUE)

  SI <- cbind(-diag(nX)*mean(subsetnew*weights), SI.beta)
  oI <- cbind(matrix(0, nrow=npar, ncol=nX),
              -t(fit$x)%*%(weights*fit$d.res)/n)
  I <- rbind(SI, oI)
  V <- (solve(I)%*%J%*%t(solve(I))*ncluster/n^2)[1:nX, 1:nX]
  vcov <- V

  out <- list(call=call, input=input, est=est, vcov=vcov)

  #---OUTPUT---

  class(out) <- "stdGee"
  return(out)
}

#' @title Summarizes GEE regression standardization fit
#'
#' @description This is a \code{summary} method for class \code{"stdGee"}.
#'
#' @param object an object of class \code{"stdGee"}.
#' @param CI.type string, indicating the type of confidence intervals. Either
#' "plain", which gives untransformed intervals, or "log", which gives
#' log-transformed intervals.
#' @param CI.level desired coverage probability of confidence intervals, on
#' decimal form.
#' @param transform a string. If set to \code{"log"}, \code{"logit"}, or
#' \code{"odds"}, the standardized mean \eqn{\theta(x)} is transformed into
#' \eqn{\psi(x)=log\{\theta(x)\}},
#' \eqn{\psi(x)=log[\theta(x)/\{1-\theta(x)\}]}, or
#' \eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively. If left unspecified,
#' \eqn{\psi(x)=\theta(x)}.
#' @param contrast a string. If set to \code{"difference"} or \code{"ratio"},
#' then \eqn{\psi(x)-\psi(x_0)} or \eqn{\psi(x) / \psi(x_0)} are constructed,
#' where \eqn{x_0} is a reference level specified by the \code{reference}
#' argument.
#' @param reference must be specified if \code{contrast} is specified.
#' @param \dots not used.
#' @author Arvid Sjolander
#' @seealso \code{\link{stdGee}}
#' @examples
#'
#' ##See documentation for stdGee
#'
#' @rdname summary
#' @export summary.stdGee
#' @export
summary.stdGee <- function(object, CI.type="plain", CI.level=0.95,
                           transform=NULL, contrast=NULL, reference=NULL, ...){

  est <- object$est
  V <- as.matrix(object$vcov)
  nX <- length(object$input$x)

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

  if(is.factor(reference))
    reference <- as.character(reference)
  est.table <- as.matrix(cbind(est, se, conf.int), nrow=length(est), ncol=4)
  dimnames(est.table) <- list(object$input$x,
                              c("Estimate", "Std. Error", paste("lower",CI.level),
                                paste("upper",CI.level)))
  out <- c(object, list(est.table=est.table,transform=transform,
                        contrast=contrast,reference=reference))

  class(out) <- "summary.stdGee"
  return(out)
}

#' @rdname print
#' @export print.stdGee
#' @export
print.stdGee <- function(x, ...){
  print(summary(x))
}

#' @title Prints summary of GEE regression standardization fit
#'
#' @description This is a \code{print} method for class \code{"summary.stdGee"}.
#'
#' @param x an object of class \code{"summary.stdGee"}.
#' @param \dots not used.
#' @author Arvid Sjolander
#' @seealso \code{\link{stdGee}}
#' @examples
#'
#' ##See documentation for stdGee
#'
#' @rdname print
#' @export print.summary.stdGee
#' @export
print.summary.stdGee <- function(x, ...){
  cat("\nFormula: ")
  print(x$input$fit$formula)
  cat("Link function:",  summary(x$input$fit)$link,  "\n")
  cat("Exposure: ", x$input$X,  "\n")
  if(!is.null(x$transform))
    cat("Transform: ", x$transform,  "\n")
  if(!is.null(x$contrast)){
    cat("Reference level: ", x$input$X, "=", x$reference,  "\n")
    cat("Contrast: ", x$contrast,  "\n")
  }
  cat("\n")
  print(x$est.table, digits=3)
}

#' @title Plots GEE regression standardization fit
#'
#' @description This is a \code{plot} method for class \code{"stdGee"}.
#'
#' @param x an object of class \code{"stdGee"}.
#' @param CI.type string, indicating the type of confidence intervals. Either
#' "plain", which gives untransformed intervals, or "log", which gives
#' log-transformed intervals.
#' @param CI.level desired coverage probability of confidence intervals, on
#' decimal form.
#' @param transform a string. If set to \code{"log"}, \code{"logit"}, or
#' \code{"odds"}, the standardized mean \eqn{\theta(x)} is transformed into
#' \eqn{\psi(x)=log\{\theta(x)\}},
#' \eqn{\psi(x)=log[\theta(x)/\{1-\theta(x)\}]}, or
#' \eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively. If left unspecified,
#' \eqn{\psi(x)=\theta(x)}.
#' @param contrast a string. If set to \code{"difference"} or \code{"ratio"},
#' then \eqn{\psi(x)-\psi(x_0)} or \eqn{\psi(x) / \psi(x_0)} are constructed,
#' where \eqn{x_0} is a reference level specified by the \code{reference}
#' argument.
#' @param reference must be specified if \code{contrast} is specified.
#' @param \dots further arguments passed on to plot.default.
#' @author Arvid Sjolander
#' @seealso \code{\link{stdGee}}
#' @examples
#'
#' ##See documentation for stdGee
#'
#' @rdname plot
#' @export plot.stdGee
#' @export
plot.stdGee <- function(x, CI.type="plain", CI.level=0.95,
                        transform=NULL, contrast=NULL, reference=NULL, ...){

  object <- x
  x <- object$input$x

  dots <- list(...)

  xlab <- object$input$X

  if(is.factor(reference))
    reference <- as.character(reference)

  if(is.null(contrast)){
    if(is.null(transform))
      ylab <- expression(mu)
    else{
      if(transform=="log")
        ylab <- expression(paste("log(", mu, ")"))
      if(transform=="logit")
        ylab <- expression(paste("logit(", mu, ")"))
      if(transform=="odds")
        ylab <- expression(paste(mu, "/(1-", mu, ")"))
    }
  }
  else{
    if(contrast=="difference"){
      if(is.null(transform))
        ylab <- c(bquote(paste(mu, "-", mu[.(reference)])), expression())
      else{
        if(transform=="log")
          ylab <- c(bquote(paste(log, "(", mu, ")-", log, "(",
                                 mu[.(reference)], ")", sep="")), expression())
        if(transform=="logit")
          ylab <- c(bquote(paste(logit, "(", mu, ")-", logit,
                                 "(", mu[.(reference)], ")", sep="")), expression())
        if(transform=="odds")
          ylab <- c(bquote(paste(mu, "/(", 1-mu, ")-",
                                 mu[.(reference)], "/(", 1-mu[.(reference)], ")", sep="")),
                    expression())
      }
    }
    if(contrast=="ratio"){
      if(is.null(transform))
        ylab <- c(bquote(paste(mu, "/", mu[.(reference)])), expression())
      else{
        if(transform=="log")
          ylab <- c(bquote(paste(log, "(", mu, ")/", log, "(",
                                 mu[.(reference)], ")", sep="")), expression())
        if(transform=="logit")
          ylab <- c(bquote(paste(logit, "(", mu, ")/", logit,
                                 "(", mu[.(reference)], ")", sep="")), expression())
        if(transform=="odds")
          ylab <- c(bquote(paste(mu, "/(", 1-mu, ")/",
                                 mu[.(reference)], "/(", 1-mu[.(reference)],
                                 ")", sep="")), expression())
      }
    }
  }

  sum.obj <- summary(object=object, CI.type=CI.type, CI.level=CI.level,
                     transform=transform, contrast=contrast, reference=reference)
  est <- sum.obj$est.table[, 1]
  lower <- sum.obj$est.table[, 3]
  upper <- sum.obj$est.table[, 4]

  ylim <- c(min(c(lower,upper)), max(c(lower,upper)))

  if(is.numeric(x) & length(x)>1){
    args <- list(x=x, y=x, xlab=xlab, ylab=ylab, ylim=ylim, type="n")
    args[names(dots)] <- dots
    do.call("plot", args=args)
    lines(x, est)
    lines(x, upper, lty=3)
    lines(x, lower, lty=3)
  }
  if(is.factor(x) | is.binary(x) | (is.numeric(x) & length(x)==1)){
    args <- list(x=1:length(x), y=1:length(x), xlab=xlab, ylab=ylab,
                 xlim=c(0, length(x)+1), ylim=ylim, type="n", xaxt="n")
    args[names(dots)] <- dots
    do.call("plot", args=args)
    points(1:length(x), est)
    points(1:length(x), upper, pch=0)
    points(1:length(x), lower, pch=0)
    for(i in 1:length(x))
      lines(x=c(i, i), y=c(lower[i], upper[i]), lty="dashed")
    mtext(text=x, side=1, at=1:length(x))
  }
}
