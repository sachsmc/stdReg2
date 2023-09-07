#' Get regression standardized estimates from a glm
#'
#' @param fit The fitted object as returned by glm.
#' @param values A named list or data.frame specifying the variables and values
#'               at which marginal means of the outcome will be estimated.
#' @param se If \code{TRUE}, calculate standard errors and the asymptotic covariance matrix.
#'
#' @returns An object of class std. This is basically a list with components estimates and covariance
#' @examples
#'
#' n <- 1000
#' Z <- rnorm(n)
#' X <- rnorm(n, mean=Z)
#' Y <- rbinom(n, 1, prob=(1+exp(X+Z))^(-1))
#' dd <- data.frame(Z, X, Y)
#' standardize_glm(formula=Y~X*Z,family="binomial", data=dd, values = list(X = 0:1))

## Tests
## Documentation

##' @importFrom stats glm
##' @export standardize_glm
standardize_glm <- function(formula.outcome, formula.exposure = NULL, data, family.outcome = gaussian, family.exposure = binomial, values, clusterid, case.control=FALSE, p.population,matched.density.cases,matched.density.controls,matching.variable) {
  n <- nrow(data)

  if (!inherits(values,c("data.frame","list"))){
    stop("values is not an object of class list or data.frame")
  }

  if(!is.data.frame(values)) {
    valuesout <- expand.grid(values)
  }
  else {
    valuesout <- values
  }
  exposure.names <- colnames(valuesout)


  ## Check that the names of values appear in the data
  xnms <- names(values)
  fitnms <- names(data)
  mnams <- match(xnms, fitnms)
  if(any(is.na(mnams))) {
    stop(paste0("variable(s) ", paste(xnms[which(is.na(mnams))], collapse = ", "),
                " not found in ", deparse1(substitute(fit.outcome)), "$data"))
  }

  if(!is.null(formula.exposure)){
    ## various checks
    if (case.control){
      warning("the setting 'case.control = TRUE' does nothing when using the doubly robust estimator")
    }
    if (!missing(clusterid)){
      warning("the option 'clusterid' does nothing when using the doubly robust estimator")
    }
    outcome.name <- as.character(formula.outcome)[2]
    exposure <- data[,exposure.names]
    if (length(exposure.names) > 1){
      stop("there has to be only when exposure with the doubly robust estimator")
    }
    ## Check that exposure is binary
    if (!(identical(family.exposure,binomial)) || !is.binary(exposure)){
      stop("the exposure has to be binary (0/1)")
    }

    fit.exposure <- tryCatch({
      stats::glm(formula=formula.exposure,data=data,family=family.exposure)
    },
    error = function(cond) {
      return(cond)
    })
    if (inherits(fit.exposure, "simpleError")){
      stop(paste0("glm for exposure failed with error: ", fit.exposure$message))
    }
    g.weights <- predict(fit.exposure, type = "response")
    weights <- exposure/g.weights+(1-exposure)/(1-g.weights)
  }
  else if (case.control){
    if (missing(p.population)){
      stop("you have to specify population prevalence when case.control = TRUE")
    }
    outcome <- as.character(formula.outcome)[2]
    outcome.data <- data[, outcome]
    if (!(identical(family.outcome,binomial)) || !is.binary(outcome.data)){
      stop("the option 'case.control = TRUE' only works with a binary outcome")
    }
    cases <- which(outcome.data==1)
    controls <- which(outcome.data==0)
    n1 <- length(cases)
    p.star <- n1/n
    n0 <- n-n1
    weights <- outcome.data*p.population/p.star+
      (1-outcome.data)*(1-p.population)/(1-p.star)*
      matched.density.controls(matching.variable)/matched.density.cases(matching.variable)
  }
  else {
    weights <- rep(1, n)
  }
  data$weights <- weights
  ## try fitting a glm model
  fit.outcome <- tryCatch({
    stats::glm(formula=formula.outcome,data=data,family=family.outcome, weights=weights)
  },
  error = function(cond) {
    return(cond)
  })
  if (inherits(fit.outcome, "simpleError")){
    stop(paste0("glm for outcome failed with error: ", fit.outcome$message))
  }

  if (is.null(formula.exposure)){
    weights <- fit.outcome$prior.weights

    ## Estimation and variance estimation
    ## In the implementation for the Score equations,
    ## the Score equations corresponding to the standardized mean come first,
    ## with the the Score equations corresponding for the parameters of the glm coming afterwards
    deriv.inv.link <- fit.outcome$family$mu.eta ## Get the derivative of the inverse link function
    predmat <- matrix(nrow = n, ncol = nrow(valuesout)) ## Contains predictions for data where the values of the data has been replaced by the i'th row of valuesout
    I.beta.means <- matrix(nrow = nrow(valuesout), ncol = length(fit.outcome$coefficients)) ## Is used for the variance calculation, where it corresponds to the derivative of the EE of theta(x) wrt. beta
    for(i in 1:nrow(valuesout)){
      ## Get the data for use with the predictions, i.e., the data from the fit but the with the variables in values replaced with the i'th row of valuesout
      data.x <- do.call("transform", c(list(data),
                                       valuesout[i,,drop = FALSE]))

      ## Save the predictions for data.x
      pred.x <- predict(object = fit.outcome, newdata = data.x)
      predmat[, i] <- fit.outcome$family$linkinv(pred.x)

      I.beta.means[i, ] <- colMeans(weights*deriv.inv.link(pred.x)*model.matrix(object=terms(fit.outcome), data=data.x)) ## Corresponds to the derivative of eta^{-1}(h(X,Z; beta)) wrt. beta, i.e., use chain rule
    }
    estimates <- colSums(weights*predmat, na.rm=TRUE)/sum(weights) ## Estimates of standardized means

    ## Implement variance estimation according to Appendix 1 of Sjolander, A. (2016)

    ## Get Score (U) and the Fisher information matrix (I) from the glm object
    sandwich.fit <- sandwich(fit=fit.outcome, data=data, weights=weights) ## NOTE: The dispersion parameter is set to 1 for the variance calculations

    ## Estimate the term that "corresponds" to var{U_{v,i}(\nu)} of Equation (5)
    EE.beta <- sandwich.fit$U
    EE.means <- weights * (predmat - matrix(rep(estimates, each=n), nrow=n,
                                            ncol=nrow(valuesout))) ## EE corresponding to the standardized means (or rather each term in estimating equation)
    EE <- cbind(EE.means, EE.beta)
    if (missing(clusterid)){
      if (case.control){
        J <- n0/n*var(EE[controls, ], na.rm=TRUE)+n1/n*var(EE[cases, ], na.rm=TRUE)
      }
      else {
        J <- var(EE, na.rm=TRUE)
      }
    }
    else {
      clusters<-data[, clusterid]
      EE <- aggr(x=EE, clusters=clusters)
      if (case.control){
        warning("case.control = TRUE may not give reasonable results for the variance with clustering")
        ## Maybe we need adjust the variance as we do above with no clustering?
      }
      J <- var(EE, na.rm=TRUE)
    }

    ## Estimate the term (I) that "corresponds" to E{\frac{\partialdU_{v,i}(\nu)}{\partial \nu}} of Equation (5)
    ## The block matrix decomposition of I can be seen in LaTeX below
    ## I=\left[
    ##   \begin{array}{ccccc|c}
    ##   -1 & 0 & 0 & \cdots & 0 &  \\
    ##   0 & -1 & 0 & \cdots & 0 &  \\
    ##   0 & 0 & -1 & \cdots & 0 & I_\beta^{(means)}  \\
    ##   \vdots & \vdots & \vdots & \ddots & \vdots &  \\
    ##   0 & 0 & 0 & \cdots & -1 & \\
    ##   \hline & & 0 & & & I^{(glm)}_\beta
    ##   \end{array}
    ## \right]

    upper.I <- cbind(-diag(nrow(valuesout))*mean(weights), I.beta.means)
    lower.I <- cbind(matrix(0, nrow=length(fit.outcome$coefficients),
                            ncol=nrow(valuesout)), sandwich.fit$I)
    I <- rbind(upper.I, lower.I)

    ## Apply Equation (5) of Sjolander, A. (2016)
    if(missing(clusterid))
      variance <- (solve(I)%*%J%*%t(solve(I))/n)[1:nrow(valuesout), 1:nrow(valuesout)]
    else
      variance <- (solve(I)%*%J%*%t(solve(I))*length(unique(clusters))/n^2)[1:nrow(valuesout), 1:nrow(valuesout)]
    fit.exposure <- NULL
  }
  else {
    ## create data
    data1 <- data0 <- data
    data1[[exposure.names]] <- 1
    data0[[exposure.names]] <- 0
    n <- nrow(data1)
    est1f <- mean(predict(fit.outcome, newdata = data1, type = "response"))
    est0f <- mean(predict(fit.outcome, newdata = data0, type = "response"))
    # refit using unweighted estimating equations
    ofitunwt <- glm(formula.outcome, family = family.outcome, data = data)
    phat <- predict(fit.exposure, type = "response")
    est1 <- predict(ofitunwt, newdata = data1, type = "response")
    est0 <- predict(ofitunwt, newdata = data0, type = "response")
    XXw <- model.matrix(fit.exposure)
    XXo <- XXo1 <- XXo0 <- model.matrix(ofitunwt)
    XXo1[,exposure.names] <- 1
    XXo0[,exposure.names] <- 0

    ## IF for the parametric models
    ifweight <- t(vcov(fit.exposure) %*% t(sandwich::estfun(fit.exposure)))
    ifout <- t(vcov(fit.outcome) %*% t(sandwich::estfun(ofitunwt)))

    eifterms1 <- (data[[exposure.names]] / phat * (data[[outcome.name]] - est1) +
                    (est1 - est1f)) / n
    eifterms0 <- ((1 - data[[exposure.names]]) / (1 - phat) * (data[[outcome.name]] - est0) +
                    (est0 - est0f)) / n ## why are we dividing by n; instead of doing it later
    hdot <- family(fit.exposure)$mu.eta(predict(fit.exposure, type = "link"))
    gdot0 <- family(fit.outcome)$mu.eta(predict(ofitunwt, newdata = data0, type = "link"))
    gdot1 <- family(fit.outcome)$mu.eta(predict(ofitunwt, newdata = data1, type = "link"))
    ## NOTE: chain rule is used for gi.dot and ri.dot
    Kterm1 <- (-1/n) * matrix(((data[[exposure.names]] * hdot) / phat^2) *
                                (data[[outcome.name]] - est1), nrow = 1, ncol = n) %*% XXw ## why multiply by XXw at the end?
    Kterm0 <- (1/n) * matrix((((1 - data[[exposure.names]]) * hdot) / (1 - phat)^2) *
                               (data[[outcome.name]] - est0), nrow = 1, ncol = n) %*% XXw #different sign
    Lterm1 <- (1/ n) * matrix(gdot1 * (1 - data[[exposure.names]]/phat),
                              nrow = 1, ncol = n) %*% XXo1 ## these matrices, why are they there?

    Lterm0 <- (1/ n) * matrix(gdot0 * ((1 - data[[exposure.names]])/(1 - phat) - 1),
                              nrow = 1, ncol = n) %*% XXo0 ## different sign
    fullif1 <- rowSums(cbind(eifterms1, (ifweight %*% t(Kterm1)), (ifout %*% t(Lterm1))))
    fullif0 <- rowSums(cbind(eifterms0, (ifweight %*% t(Kterm0)), (ifout %*% t(Lterm0))))
    covar <- sum(fullif0*fullif1)
    variance <- matrix(c(sum(fullif0^2),covar,covar,sum(fullif1^2)), nrow=2)
    exposure.names <- colnames(valuesout)
    estimates <- c(est0f,est1f)
  }
  ## Add names to asymptotic covariance matrix
  rownames(variance) <- colnames(variance) <-
    do.call("paste", c(lapply(1:ncol(valuesout), function(i) {
      paste(names(valuesout)[i], "=", round(valuesout[[i]], 2))
    }), sep = "; "))

  valuesout$se <- sqrt(diag(variance))
  valuesout$estimates <- estimates
  structure(list(estimates = valuesout, covariance = variance, fit.outcome=fit.outcome, fit.exposure=fit.exposure, exposure.names = exposure.names),
            class = "stdGLM")
}

##' @rdname summary
##' @export summary.stdGLM
##' @export
summary.stdGLM <- function(object, CI.type="plain", CI.level=0.95,
                           transform=NULL, contrast=NULL, reference=NULL, ...){
  est.old.table <- object$estimates
  est <- est.old.table$estimates
  V <- as.matrix(object$covariance)
  nX <- nrow(est.old.table)
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
    reference <- gsub(" ", "", reference)
    est.old.table$exposure <- do.call(paste, list(est.old.table[,object$exposure.names], sep=","))

    referencepos <- match(reference, est.old.table$exposure)
    if(is.na(referencepos))
      stop("reference must be a value in x")
    if(contrast=="difference"){
      dcontrast.dtransform <- diag(nX)
      dcontrast.dtransform[referencepos, ] <- -1
      dcontrast.dtransform[referencepos, referencepos] <- 0
      est <- est-est[referencepos]
    }
    else if(contrast=="ratio"){
      dcontrast.dtransform <- diag(1/est[referencepos], nrow=nX, ncol=nX)
      dcontrast.dtransform[referencepos, ] <- -est/est[referencepos]^2
      dcontrast.dtransform[referencepos, referencepos] <- 1
      est <- est/est[referencepos]
    }
    else {
      stop("contrast not supported.")
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
  if (!is.null(contrast)){
    est.table <- cbind(est.old.table$exposure, as.matrix(cbind(est, se, conf.int), nrow=length(est), ncol=4))
    est.table <- as.data.frame(est.table)
    colnames(est.table) <- c("Exposure","Estimate", "Std. Error", paste("lower",CI.level),paste("upper",CI.level))
  }
  else {
    est.table <- cbind(est.old.table[,object$exposure.names], as.matrix(cbind(est, se, conf.int), nrow=length(est), ncol=4))
    est.table <- as.data.frame(est.table)
    colnames(est.table) <- c(object$exposure.names,"Estimate", "Std. Error", paste("lower",CI.level),paste("upper",CI.level))
  }
  rownames(est.table) <- NULL
  out <- c(object, list(est.table=est.table,transform=transform,
                        contrast=contrast,reference=reference))

  class(out) <- "summary.stdGLM"
  return(out)
}

##' @rdname print
##' @export print.summary.stdGLM
##' @export
print.summary.stdGLM <- function(x, ...){
  cat("\nExposure formula: ")
  print(x$fit.exposure$formula)
  cat("Outcome link function:",  x$fit.exposure$family$link,  "\n")
  cat("Outcome formula: ")
  print(x$fit.outcome$formula)
  cat("Outcome family:",  x$fit.outcome$family$family,  "\n")
  cat("Outcome link function:",  x$fit.outcome$family$link,  "\n")
  cat("Exposure: ", paste(x$exposure.names,collapse=", "),  "\n")
  if(!is.null(x$transform))
    cat("Transform: ", x$transform,  "\n")
  if(!is.null(x$contrast)){
    cat("Reference level: ", x$input$X, "=", x$reference,  "\n")
    cat("Contrast: ", x$contrast,  "\n")
  }
  cat("\n")
  print(x$est.table, digits=3)
}

##' @rdname print
##' @export print.stdGLM
##' @export
print.stdGLM <- function(x, CI.level=0.95, ...){
  print(summary(x,CI.level=CI.level))
}

##' @rdname plot
##' @export plot.stdGLM
##' @export
plot.stdGLM <- function(x, CI.type="plain", CI.level=0.95,
                        transform=NULL, contrast=NULL, reference=NULL, ...){
  object <- x
  x <- object$estimates[,x$exposure.names]

  dots <- list(...)

  xlab <- object$exposure.names

  if (length(xlab) > 1){
    stop("cannot do plot with multiple exposures")
  }

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
  est <- sum.obj$est.table[, 2]
  lower <- sum.obj$est.table[, 4]
  upper <- sum.obj$est.table[, 5]
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

