globalVariables(".SD")
globalVariables(".")

aggr <- function(x, clusters){

  temp <- data.table(x)
  temp <- as.matrix(temp[, j=lapply(.SD, sum), by=clusters])[, -1]

}

confintall <- function(object, parm, level=0.95, fun, type="plain", ...){

  est <- do.call(what=fun, args=list(est=object$est))
  var <- delmet(fun=fun, est=object$est, vcov=object$vcov)
  ci <- CI(est=est, var=var, CI.type=type, CI.level=level)
  return(ci)
}

CI <- function(est, var, CI.type="plain", CI.level=0.95){

  se <- sqrt(var)
  qqq <- abs(qnorm((1-CI.level)/2))

  if(CI.type=="plain"){
    lower <- est-qqq*se
    upper <- est+qqq*se
  }
  if(CI.type=="log"){
    lower <- est*exp(-qqq*se/est)
    upper <- est*exp(qqq*se/est)
  }

  ci <- cbind(lower, upper)
  return(ci)

}

confint.stdCoxph <- confintall
confint.stdGee <- confintall
confint.stdGlm <- confintall
confint.stdParfrailty <- confintall

delmet <- function(fun, est, vcov){

  if(!is.list(vcov)){
    est <- matrix(est, nrow=1)
    vcov <- list(vcov)
  }
  p <- length(vcov)
  var <- vector(length=p)
  for(i in 1:p){
    gradient <- matrix(grad(func=fun, x=est[i, , drop=FALSE]), nrow=1)
    var[i] <- gradient%*%vcov[[i]]%*%t(gradient)
  }
  return(var)

}

expand <- function(x, names){
  if(is.vector(x))
    x <- x[match(names, names(x))]
  if(is.matrix(x))
    x <- x[match(names, rownames(x)), , drop=FALSE]
  return(x)
}

Hfun <- function(fit, data, fit.detail){
  if(inherits(x=fit, what="survfit")){
    #need to use summary(fit), since n.events and n.risk from fit
    #behave strange when there is left truncation
    ss <- summary(fit)
    strata <- ss$strata
    #n.strata is unweighted
    n.strata <- summary(strata)
    K <- length(n.strata)
    names.strata <- names(n.strata)
    time <- ss$time
    #n.event and n.risk are weighted
    n.event <- ss$n.event
    n.risk <- ss$n.risk
    dH <- n.event/n.risk
    H <- list()
    breaks <- c(0, cumsum(n.strata))
    for(k in 1:K){
      incl <- (breaks[k]+1):breaks[k+1]
      H[[k]] <- stepfun(time[incl], c(0, cumsum(dH[incl])))
    }
    names(H) <- names.strata
  }
  if(inherits(x=fit, what="coxph")){
    time <- fit.detail$time
    #dH is weighted
    dH <- fit.detail$hazard
    H <- stepfun(time, c(0, cumsum(dH)))
  }

  return(H)
}

is.binary <- function(v){

  if(is.numeric(v) & all(v==0 | v==1, na.rm=TRUE))
    TRUE
  else
    FALSE

}

logit <- function(x) log(x)-log(1-x)

odds <- function(x) x/(1-x)

parfrailty <- function(formula, data, clusterid, init){

  #---HELP FUNCTIONS---

  #likelihood
  like <- function(par){

    alpha <- exp(par[1])
    eta <- exp(par[2])
    phi <- exp(par[3])
    beta <- as.matrix(par[4:npar])
    B <- as.vector(X%*%beta)
    h <- delta*log(eta*t^(eta-1)/alpha^eta*exp(B))
    H <- (t/alpha)^eta*exp(B)
    Hstar <- (tstar/alpha)^eta*exp(B)
    h <- aggr(h, clusters)
    H <- aggr(H, clusters)
    Hstar <- aggr(Hstar, clusters)
    G <- d*log(phi)+cumsum(c(0, log(1/phi+j)))[d+1]
    ll <- sum(G+h+1/phi*log(1+phi*Hstar)-(1/phi+d)*
                log(1+phi*H))

    return(ll)

  }

  #cluster-specific score contributions
  scorefunc <- function(par){

    alpha <- exp(par[1])
    eta <- exp(par[2])
    phi <- exp(par[3])
    beta <- as.matrix(par[4:npar])

    #construct elements for gradient
    B <- as.vector(X%*%beta)
    h.eta <- delta*(1+eta*(log(t)-log(alpha)))
    h.beta <- X*delta
    H <- (t/alpha)^eta*exp(B)
    Hstar <- (tstar/alpha)^eta*exp(B)
    H.eta <- eta*log(t/alpha)*H
    Hstar.eta <- eta*log(tstar/alpha)*Hstar
    Hstar.eta[tstar==0] <- 0
    H.beta <- X*H
    Hstar.beta <- X*Hstar

    #aggregate elements over clusterid
    h.alpha <--d*eta
    h.eta <- aggr(h.eta, clusters)
    h.beta <- aggr(h.beta, clusters)
    H <- aggr(H, clusters)
    Hstar <- aggr(Hstar, clusters)
    H.alpha <--eta*H
    H.eta <- aggr(H.eta, clusters)
    Hstar.eta <- aggr(Hstar.eta, clusters)
    H.beta <- aggr(H.beta, clusters)
    Hstar.beta <- aggr(Hstar.beta, clusters)

    Hstar.alpha <--eta*Hstar
    K <- H/(1+phi*H)
    Kstar <- Hstar/(1+phi*Hstar)
    G.phi <- d-cumsum(c(0, 1/(1+phi*j)))[d+1]

    #first derivatives of the log-likelihood
    dl.dlogalpha <- h.alpha+Hstar.alpha/(1+phi*Hstar)-
      (1+phi*d)*H.alpha/(1+phi*H)
    dl.dlogeta <- h.eta+Hstar.eta/(1+phi*Hstar)-(1+phi*d)*
      H.eta/(1+phi*H)
    dl.dlogphi <- G.phi+1/phi*(log(1+phi*H)-
                                 log(1+phi*Hstar))+Kstar-(1+d*phi)*K
    dl.dlogbeta <- as.matrix(h.beta+Hstar.beta/(1+phi*Hstar)-
                               (1+phi*d)*H.beta/(1+phi*H))

    #score contributions
    scores <- cbind(dl.dlogalpha, dl.dlogeta, dl.dlogphi, dl.dlogbeta)
    return(scores)

  }

  #gradient
  gradientfunc <- function(par){

    return(colSums(scorefunc(par)))

  }

  #hessian
  hessianfunc <- function(par){

    alpha <- exp(par[1])
    eta <- exp(par[2])
    phi <- exp(par[3])
    beta <- as.matrix(par[4:npar])

    #construct elements for hessian
    B <- as.vector(X%*%beta)
    XX <- c(X)*X[rep(1:nrow(X), nbeta), ]
    h.eta <- delta*(1+eta*(log(t)-log(alpha)))
    H <- (t/alpha)^eta*exp(B)
    Hstar <- (tstar/alpha)^eta*exp(B)
    H.eta <- eta*log(t/alpha)*H
    Hstar.eta <- eta*log(tstar/alpha)*Hstar
    Hstar.eta[tstar==0] <- 0
    H.eta.eta <- H.eta+eta^2*(log(t/alpha))^2*H
    Hstar.eta.eta <- Hstar.eta+eta^2*(log(tstar/alpha))^2*Hstar
    Hstar.eta.eta[tstar==0] <- 0
    H.eta.beta <- eta*log(t/alpha)*(H*X)
    Hstar.eta.beta <- eta*log(tstar/alpha)*(Hstar*X)
    Hstar.eta.beta[tstar==0] <- 0
    H.beta <- cbind(H[rep(1:length(H), nbeta)]*XX, clusters)
    Hstar.beta <- cbind(Hstar[rep(1:length(H), nbeta)]*XX, clusters)

    #aggregate over clusterid
    h.eta <- aggr(h.eta, clusters)
    H <- aggr(H, clusters)
    Hstar <- aggr(Hstar, clusters)
    H.eta <- aggr(H.eta, clusters)
    Hstar.eta <- aggr(Hstar.eta, clusters)
    H.eta.eta <- aggr(H.eta.eta, clusters)
    Hstar.eta.eta <- aggr(Hstar.eta.eta, clusters)
    H.eta.beta <- aggr(H.eta.beta, clusters)
    Hstar.eta.beta <- aggr(Hstar.eta.beta, clusters)

    h.alpha.alpha <- 0
    h.alpha.eta <--d*eta
    h.eta.eta <- h.eta-d
    H.alpha <--eta*H
    Hstar.alpha <--eta*Hstar
    H.alpha.alpha <- eta^2*H
    Hstar.alpha.alpha <- eta^2*Hstar
    H.alpha.eta <--eta*(H+H.eta)
    Hstar.alpha.eta <--eta*(Hstar+Hstar.eta)

    K <- H/(1+phi*H)
    Kstar <- Hstar/(1+phi*Hstar)
    G.phi.phi <- cumsum(c(0,phi*j/(1+phi*j)^2))[d+1]

    #derivative of gradient wrt logalpha wrt all parameters except beta
    dl.dlogalpha.dlogalpha <- sum(h.alpha.alpha+Hstar.alpha.alpha/
                                    (1+phi*Hstar)-phi*(Hstar.alpha/(1+phi*Hstar))^2 -
                                    (1+phi*d)*(H.alpha.alpha/(1+phi*H) -
                                                 phi*(H.alpha/(1+phi*H))^2))
    dl.dlogalpha.dlogeta <- sum(h.alpha.eta+Hstar.alpha.eta/
                                  (1+phi*Hstar)-phi*Hstar.alpha*Hstar.eta/
                                  (1+phi*Hstar)^2-(1+phi*d)*(H.alpha.eta/(1+phi*H) -
                                                               phi*H.alpha*H.eta/(1+phi*H)^2))
    dl.dlogalpha.dlogphi <- sum(phi*(-Hstar.alpha*Hstar/
                                       (1+phi*Hstar)^2+H.alpha*H/(1+phi*H)^2-
                                       d*(H.alpha/(1+phi*H)-phi*H.alpha*H/(1+phi*H)^2)))
    ddl.dlogalpha <- cbind(dl.dlogalpha.dlogalpha, dl.dlogalpha.dlogeta,
                           dl.dlogalpha.dlogphi)

    #derivative of gradient wrt logeta wrt all parameters except beta
    dl.dlogeta.dlogeta <- sum(h.eta.eta+Hstar.eta.eta/(1+phi*Hstar) -
                                phi*(Hstar.eta/(1+phi*Hstar))^2-(1+phi*d) *
                                (H.eta.eta/(1+phi*H)-phi*(H.eta/(1+phi*H))^2))
    dl.dlogeta.dlogphi <- sum(phi*(-Hstar.eta*Hstar/
                                     (1+phi*Hstar)^2+H.eta*H/(1+phi*H)^2-
                                     d*(H.eta/(1+phi*H)-phi*H.eta*H/(1+phi*H)^2)))
    ddl.dlogeta <- cbind(dl.dlogalpha.dlogeta, dl.dlogeta.dlogeta,
                         dl.dlogeta.dlogphi)

    #derivative of gradient wrt logphi wrt all parameters except beta
    dl.dlogphi.dlogphi <- sum(G.phi.phi+1/phi*
                                (log(1+phi*Hstar)-log(1+phi*H))+ K-Kstar+phi
                              *(K^2-Kstar^2)+d*phi*K*(phi*K-1))
    ddl.dlogphi <- cbind(dl.dlogalpha.dlogphi, dl.dlogeta.dlogphi,
                         dl.dlogphi.dlogphi)

    #derivatives of gradients wrt (logalpha, logeta, logphi) wrt beta
    H <- (t/alpha)^eta*exp(B)
    Hstar <- (tstar/alpha)^eta*exp(B)
    XX <- c(X)*X[rep(1:nrow(X), nbeta), ]
    nbeta_rep <- rep(1:nbeta, each=nrow(X))
    H.beta <- as.matrix(aggr(H*X, clusters))
    H.beta2 <- H.beta[rep(1:nrow(H.beta), nbeta), ]*c(H.beta)
    Hstar.beta <- as.matrix(aggr(Hstar*X, clusters))
    Hstar.beta2 <- Hstar.beta[rep(1:nrow(Hstar.beta), nbeta), ]*c(Hstar.beta)
    Hstar.beta.beta <- data.table(nbeta_rep, clusters, Hstar*XX)
    H.beta.beta <- data.table(nbeta_rep, clusters, H*XX)
    H <- aggr(H, clusters)
    Hstar <- aggr(Hstar, clusters)
    Hstar2 <- phi*Hstar.beta2/(1+phi*Hstar)^2
    H2 <- phi*(1+d*phi)*H.beta2/(1+phi*H)^2
    Hstar.beta.beta <- data.table(clusters, nbeta_rep, Hstar.beta.beta)
    Hstar.beta.beta <- as.matrix(Hstar.beta.beta[, j=lapply(.SD, sum),
                                                 by=.(nbeta_rep, clusters)])[, -1:-2, drop=FALSE]
    H.beta.beta <- data.table(clusters, nbeta_rep, H.beta.beta)
    H.beta.beta <- as.matrix(H.beta.beta[, j=lapply(.SD, sum),
                                         by=.(nbeta_rep, clusters)])[, -1:-2, drop=FALSE]
    Hstar1 <- Hstar.beta.beta/(1+phi*Hstar)
    H1 <- (1+d*phi)*H.beta.beta/(1+phi*H)
    H.alpha.beta <--eta*H.beta
    Hstar.alpha.beta <--eta*Hstar.beta

    dl.dlogalpha.dlogbeta <- colSums(as.matrix(Hstar.alpha.beta/
                                                 (1+phi*Hstar)-phi*Hstar.alpha*Hstar.beta/
                                                 (1+phi*Hstar)^2-(1+phi*d)*(H.alpha.beta/
                                                                              (1+phi*H)-phi*H.alpha*H.beta/(1+phi*H)^2)))
    ddl.dlogalpha <- cbind(ddl.dlogalpha, t(dl.dlogalpha.dlogbeta))


    dl.dlogeta.dlogbeta <- t(colSums(as.matrix(Hstar.eta.beta/
                                                 (1+phi*Hstar)-phi*Hstar.eta*Hstar.beta/
                                                 (1+phi*Hstar)^2-(1+phi*d)*(H.eta.beta/(1+phi*H) -
                                                                              phi*H.eta*H.beta/(1+phi*H)^2))))
    ddl.dlogeta <- cbind(ddl.dlogeta, dl.dlogeta.dlogbeta)


    dl.dlogphi.dlogbeta <- t(colSums(as.matrix(phi*
                                                 (-Hstar.beta*Hstar/(1+phi*Hstar)^2+H.beta*H/
                                                    (1+phi*H)^2-d*(H.beta/(1+phi*H)-phi*H.beta*H/
                                                                     (1+phi*H)^2)))))
    ddl.dlogphi <- cbind(ddl.dlogphi, dl.dlogphi.dlogbeta)

    #derivative of gradient wrt to beta wrt to all parameters
    dl.dlogbeta.dlogbeta <- (Hstar.beta.beta/(1+phi*Hstar) -
                               phi*Hstar.beta2/(1+phi*Hstar)^2-(1+phi*d) *
                               (H.beta.beta/(1+phi*H)-phi*H.beta2/(1+phi*H)^2))
    nbeta_rep2 <- rep(1:nbeta, each=length(H))
    dl.dlogbeta.dlogbeta <- data.table(nbeta_rep2, dl.dlogbeta.dlogbeta)
    dl.dlogbeta.dlogbeta <- as.matrix(dl.dlogbeta.dlogbeta[,
                                                           j=lapply(.SD, sum), by=.(nbeta_rep2)])[, -1]
    ddl.dlogbeta <- rbind(dl.dlogalpha.dlogbeta, dl.dlogeta.dlogbeta,
                          dl.dlogphi.dlogbeta, dl.dlogbeta.dlogbeta)

    hessian <- cbind(t(ddl.dlogalpha), t(ddl.dlogeta), t(ddl.dlogphi),
                     ddl.dlogbeta)
    return(hessian)

  }

  #---PREPARATION---
  call <- match.call()

  #delete rows with missing on variables in the model
  data.temp <- data
  m <- model.matrix(object=formula, data=data.temp)
  data.temp <- data.temp[match(rownames(m), rownames(data.temp)), ]
  X <- model.matrix(formula, data=data.temp)[, -1, drop=FALSE]
  clusters <- data.temp[, clusterid]
  n <- nrow(X)
  ncluster <- length(unique(clusters))
  nbeta <- ncol(X)
  npar <- 3+nbeta
  if(missing(init)) init <- c(rep(0, npar))

  #extract start variable, end variable and event variable
  Y <- model.extract(frame=model.frame(formula=formula, data=data.temp),
                     component="response")
  if (ncol(Y)==2) {
    tstar <- rep(0, nrow(data.temp))
    t <- Y[, 1]
    delta <- Y[, 2]
  }
  if (ncol(Y)==3) {
    tstar <- Y[, 1]
    t <- Y[, 2]
    delta <- Y[, 3]
  }

  #number of uncensored in each cluster
  d <- aggr(delta, clusters)
  D <- max(d)
  j <- 0:(D-1)

  #---MODEL FITTING

  #maximize log likelihood
  fit <- optim(par=init, fn=like, gr=gradientfunc, method="BFGS", hessian=FALSE,
               control=list(fnscale=-1))
  est <- fit$par
  names(est) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
                  colnames(X))

  #calculate score contributions
  score <- scorefunc(par=est)
  colnames(score) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
                       colnames(X))

  #calculate hessian
  hessian <- hessianfunc(par=est)
  rownames(hessian) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
                         colnames(X))
  colnames(hessian) <- c("log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
                         colnames(X))

  output <- c(list(formula=formula, data=data, clusterid=clusterid,
                   ncluster=ncluster, n=n, X=X, fit=fit, est=est,
                   score=score, vcov=-solve(hessian), call=call))

  class(output) <- "parfrailty"
  return(output)
}

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

plot.stdParfrailty <- plot.stdCoxph

plot.stdGlm <- function(x, CI.type="plain", CI.level=0.95,
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

plot.stdGee <- plot.stdGlm

print.parfrailty <- function(x, ...){

  cat("Call:", "\n")
  print.default(x$call)
  cat("\nEstimated parameters in the Gamma-Weibull frailty model", "\n")
  cat("\n")
  table.est <- as.matrix(cbind(x$est, sqrt(diag(x$vcov))))
  rownames(table.est) <- names(x$est)
  colnames(table.est) <- c("coef", "se(coef)")
  print.default(table.est)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of clusters:", x$ncluster)
  cat("\n")

}

print.summary.parfrailty <- function(x, digits=max(3L, getOption("digits")-3L),
                                     ...){
  ## Function call
  cat("Call:", "\n")
  print.default(x$call)
  level <- x$CI.level*100
  CI.text <- paste0(as.character(level),"%")
  cat("\nEstimated parameters in the Gamma-Weibull frailty model", "\n")
  cat("\n")
  table.est <- cbind(x$est, exp(x$est), x$se, x$zvalue, x$pvalue)
  rownames(table.est) <- names(x$est)
  colnames(table.est) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  #print.default(table.est, digits=3)
  printCoefmat(table.est, digits=3)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of clusters:", x$ncluster, "\n")
  cat("\n")

}

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

print.summary.stdParfrailty <- print.summary.stdCoxph

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

print.summary.stdGlm <- function(x, ...){

  cat("\nFormula: ")
  print(x$input$fit$formula)
  cat("Family:",  x$input$fit$family$family,  "\n")
  cat("Link function:",  x$input$fit$family$link,  "\n")
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

sandwich <- function(fit, data, weights, t, fit.detail){

  n <- nrow(data)
  if(missing(weights))
    weights <- rep(1, n)

  if(inherits(x=fit, what="glm")){

    #---meat---

    m <- expand(model.matrix(fit), rownames(data))
    res <- expand(residuals(fit, type="response"), rownames(data))
    U <- weights*m*res
    U[is.na(U)] <- 0

    #---bread---

    #summary(fit)$cov.unscaled is weighted
    I <- -solve(summary(fit)$cov.unscaled)/n

  }
  if(inherits(x=fit, what="ah")){

    #---meat---

    #residuals are weighted
    res <- predict(object=fit, type="residuals")
    rownames(res) <- fit$incl
    colnames(res) <- names(fit$coefficients)
    res <- expand(res, rownames(data))
    U <- res

  }
  if(inherits(x=fit, what="coxph")){

    #---meat for regression coefficients---

    #score residuals are unweighted, but computed with estimated coefficients
    #from weighted model
    res <- residuals(fit, type="score")
    res <- expand(res, rownames(data))
    Ucoef <- weights*res

    #---bread for regression coefficients---

    #vcov(fit) is weighted
    Icoef <- -solve(vcov(fit))/n

    if(missing(t)){

      U <- Ucoef
      colnames(U) <- names(fit$coef)
      U[is.na(U)] <- 0
      I <- Icoef

    }else{

      nt <- length(t)

      #---meat and bread for baseline hazard---

      #check if left truncation
      varsLHS <- all.vars(fit$formula[[2]])
      t2 <- varsLHS[length(varsLHS)-1]
      if(length(varsLHS)==3)
        t1 <- varsLHS[1]
      else
        t1 <- NULL
      time <- fit.detail$time
      #nevent is unweighted
      nevent <- fit.detail$nevent
      #hazard is weighted
      dH <- fit.detail$hazard
      names(dH) <- time
      #varhaz is dH/mean(exp(x*b)) where mean is in risk set
      #varhaz is weighted
      dHvar <- fit.detail$varhaz
      Hvar <- stepfun(time, c(0, cumsum(dHvar)))
      p <- predict(object=fit, type="risk")
      m <- model.matrix(fit)
      names(p) <- rownames(m)
      p <- expand(p, rownames(data))
      ncoef <- length(fit$coef)
      #fit.detail$means is weighted, but not relative to the mean covariate
      #in the sample, like all other outputs from coxph.detail,
      #subtracting fit$means fixes this
      means <- as.matrix(fit.detail$means)
      means <- means-matrix(fit$means, nrow=nrow(means), ncol=ncol(means),
                            byrow=TRUE)
      means <- expand(means, data[, t2])
      UH <- matrix(nrow=n, ncol=nt)
      IH <- matrix(nrow=nt, ncol=ncoef)
      for(j in 1:nt){
        #dividing with nevent accounts for ties,
        #but this assumes that H is computed with Breslow method for ties,
        #and that weights are equal within ties.
        tmp1 <- n*expand(dH/nevent*(time<=t[j]), data[, t2])
        tmp1[is.na(tmp1)] <- 0
        tmp2 <- n*Hvar(pmin(t[j], data[, t2]))*p
        if(!is.null(t1))
          tmp2 <- tmp2-n*Hvar(data[, t1])*(data[, t1]<t[j])*p
        UH[, j] <- tmp1-tmp2
        dH.dbeta <- means*tmp1
        dH.dbeta[is.na(dH.dbeta)] <- 0
        IH[j, ] <- -colMeans(dH.dbeta)
      }
      U <- cbind(Ucoef, UH)
      colnames(U) <- c(names(fit$coef), paste0("H", t))
      U[is.na(U)] <- 0
      I <- rbind(cbind(Icoef, matrix(0, nrow=ncoef, ncol=length(t))),
                 cbind(IH, -diag(length(t))))
    }

  }
  if(inherits(x=fit, what="survfit")){

    #---meat---

    #check if left truncation
    #survfit object has no formula element, so need to get it from call,
    #need to use eval, since the fit$call$formula will be literary what the user
    #gave as argument, e.g. if formula=f, then fit$call$formula is f, not the
    #formula contained in f
    varsLHS <- all.vars(eval(fit$call$formula)[[2]])
    t2 <- varsLHS[length(varsLHS)-1]
    if(length(varsLHS)==3)
      t1 <- varsLHS[1]
    else
      t1 <- NULL
    #need to use summary(fit), since n.events and n.risk from fit
    #behave strange when there is left truncation
    ss <- summary(fit)
    strata <- ss$strata
    #n.strata is unweighted
    n.strata <- summary(strata)
    K <- length(n.strata)
    names.strata <- names(n.strata)
    time <- ss$time
    #n.event and n.risk are weighted
    n.event <- ss$n.event
    n.risk <- ss$n.risk
    dH <- n.event/n.risk
    names(dH) <- paste(time, strata)
    dHvar <- dH/n.risk
    #survfit object has no formula element, so need to get it from call,
    #need to use eval, since the fit$call$formula will be literary what the user
    #gave as argument, e.g. if formula=f, then fit$call$formula is f, not the
    #formula contained in f
    vars <- all.vars(eval(fit$call$formula)[[3]])
    #note: strata is a function in the survival package
    strata.all <- strata(data[, vars, drop=FALSE])
    tmp1 <- matrix(nrow=n, ncol=K)
    U <- matrix(nrow=n, ncol=K)
    colnames(U) <- names.strata
    breaks <- c(0, cumsum(n.strata))
    for(k in 1:K){
      incl <- (breaks[k]+1):breaks[k+1]
      Hvar <- stepfun(time[incl], c(0, cumsum(dHvar[incl])))
      #dividing with nevent[incl] account for ties,
      #but this assumes that H is computed with Breslow method for ties,
      #and that weights are equal within ties.
      #multiplying with weights corrects for nevent being weighted;
      #here we just want to divide with the actual number of events to account
      #for ties, not the weighted number of events
      tmp1.time <- n*dH[incl]/n.event[incl]*(time[incl]<=t)
      tmp1[, k] <- tmp1.time[match(paste(data[, t2], strata.all),
                                   names(tmp1.time))]*weights
      tmp1[is.na(tmp1[, k]), k] <- 0
      sk <- names.strata[k]
      incl <- which(strata.all==sk)
      tmp2 <- n*Hvar(pmin(t, data[incl, t2]))
      if(!is.null(t1))
        tmp2 <- tmp2-n*Hvar(data[incl, t1])*(data[incl, t1]<t)
      U[incl, k] <- tmp1[incl, k]-tmp2
    }

    #---bread---

    I <- diag(-1, K)
    rownames(I) <- names.strata
    colnames(I) <- names.strata

  }

  U[is.na(U)] <- 0
  return(list(I=I, U=U))

}

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

stdGlm <- function(fit, data, X, x, clusterid, case.control=FALSE, subsetnew){

  call <- match.call()

  #---PREPARATION---

  formula <- fit$formula
  weights <- fit$prior.weights
  npar <- length(fit$coef)

  if(case.control & all(fit$prior.weights==1))
    warning("stdGlm only gives consistent estimates with case-control data,
      if weights are used that corrects for the biased sampling scheme")

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

  #---ESTIMATES OF MEANS AT VALUES SPECIFIED BY x ---

  pred <- matrix(nrow=n, ncol=nX)
  g <- family(fit)$mu.eta
  SI.beta <- matrix(nrow=nX, ncol=npar)
  for(i in 1:nX){
    data.x <- data
    if(!is.na(x[i]))
      data.x[, X] <- x[i]
    pred[, i] <- predict(object=fit, newdata=data.x, type="response")
    dmu.deta <- g(predict(object=fit, newdata=data.x))
    #Need object=terms(fit) here. If formula contains splines,
    #then neither object=formula nor object=fit will not work when no variation
    #in exposure, since model.matrix needs to retrieve Boundary.knots
    #from terms(fit).
    deta.dbeta <- model.matrix(object=terms(fit), data=data.x)
    dmu.dbeta <- dmu.deta*deta.dbeta
    SI.beta[i, ] <- colMeans(subsetnew*weights*dmu.dbeta)
  }
  est <- colSums(subsetnew*weights*pred, na.rm=TRUE)/
    sum(subsetnew*weights)

  #---VARIANCE OF MEANS AT VALUES SPECIFIED BY x---
  sandwich.fit <- sandwich(fit=fit, data=data, weights=weights)
  ores <- sandwich.fit$U
  mres <- subsetnew*weights*(pred-matrix(rep(est, each=n), nrow=n,
                                         ncol=nX))
  res <- cbind(mres, ores)

  if(case.control & missing(clusterid)){
    outcome <- as.character(formula)[2]
    controls <- which(data[, outcome]==0)
    cases <- which(data[, outcome]==1)
    n0 <- length(controls)
    n1 <- length(cases)
    J <- n0/n*var(res[controls, ], na.rm=TRUE)+n1/n*var(res[cases, ],
                                                        na.rm=TRUE)
  }
  else{
    if(!missing(clusterid))
      res <- aggr(x=res, clusters=data[, clusterid])
    J <- var(res, na.rm=TRUE)
  }

  SI <- cbind(-diag(nX)*mean(subsetnew*weights), SI.beta)
  oI <- cbind(matrix(0, nrow=npar, ncol=nX), sandwich.fit$I)
  I <- rbind(SI, oI)
  if(missing(clusterid))
    V <- (solve(I)%*%J%*%t(solve(I))/n)[1:nX, 1:nX]
  else
    V <- (solve(I)%*%J%*%t(solve(I))*ncluster/n^2)[1:nX, 1:nX]
  vcov <- V

  out <- list(call=call, input=input, est=est, vcov=vcov)

  #---OUTPUT---

  class(out) <- "stdGlm"
  return(out)

}

stdParfrailty <- function(fit, data, X, x, t, clusterid, subsetnew){

  call <- match.call()

  #---PREPARATION---

  formula <- fit$formula
  npar <- length(fit$est)

  #delete rows that did not contribute to the model fit,
  #e.g. missing data or in subset for fit
  #note: object=fit does not work when fit is parfrailty object
  m <- model.matrix(object=formula, data=data)
  data <- data[match(rownames(m), rownames(data)), ]
  n <- nrow(data)

  #Make new subset if supplied.
  subsetnew <-
    if(missing(subsetnew))
      rep(1, n)
  else
    as.numeric(eval(substitute(subsetnew), data, parent.frame()))
  input <- as.list(environment())

  #extract end variable and event variable
  Y <- model.extract(frame=model.frame(formula=formula, data=data),
                     component="response")
  if(ncol(Y)==2){
    end <- Y[, 1]
    event <- Y[, 2]
  }
  if(ncol(Y)==3){
    end <- Y[, 2]
    event <- Y[, 3]
  }

  clusters <- data[, clusterid]
  ncluster <- length(unique(clusters))

  #assign values to x and reference if not supplied
  #make sure x is a factor if data[, X] is a factor,
  #with the same levels as data[, X]
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

  #assign value to t if missing
  if(missing(t))
    t <- end[event==1]
  t <- sort(t)
  input$t <- t
  nt <- length(t)

  #preparation
  est <- matrix(nrow=nt, ncol=nX)
  vcov <- vector(mode="list", length=nt)
  logalpha <- fit$est[1]
  alpha <- exp(logalpha)
  logeta <- fit$est[2]
  eta <- exp(logeta)
  logphi <- fit$est[3]
  phi <- exp(logphi)
  beta <- fit$est[(3+1):npar]

  #---LOOP OVER nt

  for(j in 1:nt){

    if(t[j]==0){
      est[j, ] <- 1
      vcov[[j]] <- matrix(0, nrow=nX, ncol=nX)
    }
    else{

      H0t <- (t[j]/alpha)^eta

      #---ESTIMATES OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---

      si <- matrix(nrow=n, ncol=nX)
      SI.logalpha <- vector(length=nX)
      SI.logeta <- vector(length=nX)
      SI.logphi <- vector(length=nX)
      SI.beta <- matrix(nrow=nX, ncol=npar-3)
      for(i in 1:nX){
        data.x <- data
        if(!is.na(x[i]))
          data.x[, X] <- x[i]
        m <- model.matrix(object=formula, data=data.x)[, -1, drop=FALSE]
        predX <- colSums(beta*t(m))
        temp <- 1+phi*H0t*exp(predX)
        si[, i] <- temp^(-1/phi)
        SI.logalpha[i] <- mean(subsetnew*H0t*eta*exp(predX)/temp^(1/phi+1))*n/
          ncluster
        SI.logeta[i] <- mean(subsetnew*(-H0t)*exp(predX)*log(t[j]/alpha)*eta/
                               temp^(1/phi+1))*n/ncluster
        SI.logphi[i] <- mean(subsetnew*log(temp)/(phi*temp^(1/phi))-
                               H0t*exp(predX)/temp^(1/phi+1))*n/ncluster
        SI.beta[i, ] <- colMeans(subsetnew*(-H0t)*exp(predX)*m/temp^(1/phi+1))*
          n/ncluster
      }
      est[j, ] <- colSums(subsetnew*si, na.rm=TRUE)/sum(subsetnew)

      #---VARIANCE OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---

      sres <- subsetnew*(si-matrix(rep(est[j, ],each=n), nrow=n, ncol=nX))
      sres <- aggr(sres, clusters)
      coefres <- fit$score
      res <- cbind(sres, coefres)

      J <- var(res, na.rm=TRUE)

      #Note: the term n/ncluster is because SI.logalpha, SI.logeta, SI.logphi,
      #and SI.beta are clustered, which they are not in stdCoxph
      SI <- cbind(-diag(nX)*mean(subsetnew)*n/ncluster, SI.logalpha, SI.logeta,
                  SI.logphi, SI.beta)

      betaI <- cbind(matrix(0, nrow=npar, ncol=nX), -solve(fit$vcov)/ncluster)
      I <- rbind(SI, betaI)
      V <- (solve(I)%*%J%*%t(solve(I))/ncluster)[1:nX, 1:nX]
      vcov[[j]] <- V

    }

  }

  out <- list(call=call, input=input, est=est, vcov=vcov)


  #---OUTPUT---

  class(out) <- "stdParfrailty"
  return(out)

}

summary.parfrailty <- function(object, CI.type="plain", CI.level=0.95,
                               digits=max(3L, getOption("digits")-3L), ...){
  if(missing(CI.level)) CI.level <- 0.95
  if(missing(CI.type)) CI.type <- "plain"

  ### Inference
  var <- diag(object$vcov)
  se <- sqrt(var)
  zvalue <- object$est/se
  pvalue <- 2*pnorm(-abs(zvalue))
  confidence.interval <- CI(est=object$est, var=var,
                            CI.type=CI.type, CI.level=CI.level)
  colnames(confidence.interval) <- c("Lower limit", "Upper limit")

  ans <- list(est=object$est, se=se, zvalue=zvalue,
              pvalue=pvalue, score=object$score, X=object$X, vcov=object$vcov,
              call=object$call, formula=object$formula, modelmatrix=object$X,
              data=object$data, clusterid=object$clusterid,
              ncluster=object$ncluster, n=object$n, CI.type=CI.type, CI.level=CI.level,
              confidence.interval=confidence.interval)

  class(ans) <- "summary.parfrailty"
  return(ans)

}

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

summary.stdParfrailty <- summary.stdCoxph

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

summary.stdGlm <- function(object, CI.type="plain", CI.level=0.95,
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

  class(out) <- "summary.stdGlm"
  return(out)

}











