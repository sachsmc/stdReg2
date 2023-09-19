#' @title Fits shared frailty gamma-Weibull models
#'
#' @description \code{parfrailty} fits shared frailty gamma-Weibull models. It is
#' specifically designed to work with the function \code{stdParfrailty}, which
#' performs regression standardization in shared frailty gamma-Weibull models.
#'
#' @details \code{parfrailty} fits the shared frailty gamma-Weibull model
#' \deqn{\lambda(t_{ij}|C_{ij})=\lambda(t_{ij};\alpha,\eta)U_iexp\{h(C_{ij};\beta)\},}
#' where \eqn{t_{ij}} and \eqn{C_{ij}} are the survival time and covariate
#' vector for subject \eqn{j} in cluster \eqn{i}, respectively.
#' \eqn{\lambda(t;\alpha,\eta)} is the Weibull baseline hazard function
#' \deqn{\eta t^{\eta-1}\alpha^{-\eta},} where \eqn{\eta} is the shape
#' parameter and \eqn{\alpha} is the scale parameter. \eqn{U_i} is the
#' unobserved frailty term for cluster \eqn{i}, which is assumed to have a
#' gamma distribution with scale = 1/shape = \eqn{\phi}. \eqn{h(X;\beta)} is
#' the regression function as specified by the \code{formula} argument,
#' parametrized by a vector \eqn{\beta}. The ML estimates
#' \eqn{\{log(\hat{\alpha}),log(\hat{\eta}),log(\hat{\phi}),\hat{\beta}\}} are
#' obtained by maximizing the marginal (over \eqn{U}) likelihood.
#'
#' @param formula an object of class "\code{formula}", on the same format as
#' accepted by the \code{coxph} function in the \pkg{survival} package.
#' @param data a data frame containing the variables in the model.
#' @param clusterid an string containing the name of a cluster identification
#' variable.
#' @param init an optional vector of initial values for the model parameters.
#' @return An object of class \code{"parfrailty"} is a list containing:
#' \item{est}{ the ML estimates \eqn{\{log(\hat{\alpha}),log(\hat{\eta}),
#' log(\hat{\phi}),\hat{\beta}\}}. } \item{vcov}{ the variance-covariance
#' vector of the ML estimates. } \item{score}{ a matrix containing the
#' cluster-specific contributions to the ML score equations.  }
#' @note If left truncation is present, it is assumed that it is strong left
#' truncation.  This means that, even if the truncation time may be
#' subject-specific, the whole cluster is unobserved if at least one subject in
#' the cluster dies before his/her truncation time. If all subjects in the
#' cluster survive beyond their subject-specific truncation times, then the
#' whole cluster is observed (Van den Berg and Drepper, 2016).
#' @author Arvid Sjolander and Elisabeth Dahlqwist.
#' @references Dahlqwist E., Pawitan Y., Sjolander A. (2019). Regression
#' standardization and attributable fraction estimation with between-within
#' frailty models for clustered survival data. \emph{Statistical Methods in
#' Medical Research} \bold{28}(2), 462-485.
#'
#' Van den Berg G.J., Drepper B. (2016). Inference for shared frailty survival
#' models with left-truncated data. \emph{Econometric Reviews}, 35(6),
#' 1075-1098.
#' @examples
#'
#'
#' require(survival)
#'
#' # simulate data
#' set.seed(5)
#' n <- 200
#' m <- 3
#' alpha <- 1.5
#' eta <- 1
#' phi <- 0.5
#' beta <- 1
#' id <- rep(1:n, each = m)
#' U <- rep(rgamma(n, shape = 1 / phi, scale = phi), each = m)
#' X <- rnorm(n * m)
#' # reparametrize scale as in rweibull function
#' weibull.scale <- alpha / (U * exp(beta * X))^(1 / eta)
#' T <- rweibull(n * m, shape = eta, scale = weibull.scale)
#'
#' # right censoring
#' C <- runif(n * m, 0, 10)
#' D <- as.numeric(T < C)
#' T <- pmin(T, C)
#'
#' # strong left-truncation
#' L <- runif(n * m, 0, 2)
#' incl <- T > L
#' incl <- ave(x = incl, id, FUN = sum) == m
#' dd <- data.frame(L, T, D, X, id)
#' dd <- dd[incl, ]
#'
#' fit <- parfrailty(formula = Surv(L, T, D) ~ X, data = dd, clusterid = "id")
#' print(summary(fit))
#'
#' @export parfrailty
parfrailty <- function(formula, data, clusterid, init) {
  #---HELPER FUNCTIONS---

  # likelihood
  like <- function(par) {
    alpha <- exp(par[1])
    eta <- exp(par[2])
    phi <- exp(par[3])
    beta <- as.matrix(par[4:npar])
    B <- as.vector(X %*% beta)
    h <- delta * log(eta * t^(eta - 1) / alpha^eta * exp(B))
    H <- (t / alpha)^eta * exp(B)
    Hstar <- (tstar / alpha)^eta * exp(B)
    h <- aggr(h, clusters)
    H <- aggr(H, clusters)
    Hstar <- aggr(Hstar, clusters)
    G <- d * log(phi) + cumsum(c(0, log(1 / phi + j)))[d + 1]
    ll <- sum(G + h + 1 / phi * log(1 + phi * Hstar) - (1 / phi + d) *
      log(1 + phi * H))

    return(ll)
  }

  # cluster-specific score contributions
  scorefunc <- function(par) {
    alpha <- exp(par[1])
    eta <- exp(par[2])
    phi <- exp(par[3])
    beta <- as.matrix(par[4:npar])

    # construct elements for gradient
    B <- as.vector(X %*% beta)
    h.eta <- delta * (1 + eta * (log(t) - log(alpha)))
    h.beta <- X * delta
    H <- (t / alpha)^eta * exp(B)
    Hstar <- (tstar / alpha)^eta * exp(B)
    H.eta <- eta * log(t / alpha) * H
    Hstar.eta <- eta * log(tstar / alpha) * Hstar
    Hstar.eta[tstar == 0] <- 0
    H.beta <- X * H
    Hstar.beta <- X * Hstar

    # aggregate elements over clusterid
    h.alpha <- -d * eta
    h.eta <- aggr(h.eta, clusters)
    h.beta <- aggr(h.beta, clusters)
    H <- aggr(H, clusters)
    Hstar <- aggr(Hstar, clusters)
    H.alpha <- -eta * H
    H.eta <- aggr(H.eta, clusters)
    Hstar.eta <- aggr(Hstar.eta, clusters)
    H.beta <- aggr(H.beta, clusters)
    Hstar.beta <- aggr(Hstar.beta, clusters)

    Hstar.alpha <- -eta * Hstar
    K <- H / (1 + phi * H)
    Kstar <- Hstar / (1 + phi * Hstar)
    G.phi <- d - cumsum(c(0, 1 / (1 + phi * j)))[d + 1]

    # first derivatives of the log-likelihood
    dl.dlogalpha <- h.alpha + Hstar.alpha / (1 + phi * Hstar) -
      (1 + phi * d) * H.alpha / (1 + phi * H)
    dl.dlogeta <- h.eta + Hstar.eta / (1 + phi * Hstar) - (1 + phi * d) *
      H.eta / (1 + phi * H)
    dl.dlogphi <- G.phi + 1 / phi * (log(1 + phi * H) -
      log(1 + phi * Hstar)) + Kstar - (1 + d * phi) * K
    dl.dlogbeta <- as.matrix(h.beta + Hstar.beta / (1 + phi * Hstar) -
      (1 + phi * d) * H.beta / (1 + phi * H))

    # score contributions
    scores <- cbind(dl.dlogalpha, dl.dlogeta, dl.dlogphi, dl.dlogbeta)
    return(scores)
  }

  # gradient
  gradientfunc <- function(par) {
    return(colSums(scorefunc(par)))
  }

  # hessian
  hessianfunc <- function(par) {
    alpha <- exp(par[1])
    eta <- exp(par[2])
    phi <- exp(par[3])
    beta <- as.matrix(par[4:npar])

    # construct elements for hessian
    B <- as.vector(X %*% beta)
    XX <- c(X) * X[rep(1:nrow(X), nbeta), ]
    h.eta <- delta * (1 + eta * (log(t) - log(alpha)))
    H <- (t / alpha)^eta * exp(B)
    Hstar <- (tstar / alpha)^eta * exp(B)
    H.eta <- eta * log(t / alpha) * H
    Hstar.eta <- eta * log(tstar / alpha) * Hstar
    Hstar.eta[tstar == 0] <- 0
    H.eta.eta <- H.eta + eta^2 * (log(t / alpha))^2 * H
    Hstar.eta.eta <- Hstar.eta + eta^2 * (log(tstar / alpha))^2 * Hstar
    Hstar.eta.eta[tstar == 0] <- 0
    H.eta.beta <- eta * log(t / alpha) * (H * X)
    Hstar.eta.beta <- eta * log(tstar / alpha) * (Hstar * X)
    Hstar.eta.beta[tstar == 0] <- 0
    H.beta <- cbind(H[rep(seq_len(length(H)), nbeta)] * XX, clusters)
    Hstar.beta <- cbind(Hstar[rep(seq_len(length(H)), nbeta)] * XX, clusters)

    # aggregate over clusterid
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
    h.alpha.eta <- -d * eta
    h.eta.eta <- h.eta - d
    H.alpha <- -eta * H
    Hstar.alpha <- -eta * Hstar
    H.alpha.alpha <- eta^2 * H
    Hstar.alpha.alpha <- eta^2 * Hstar
    H.alpha.eta <- -eta * (H + H.eta)
    Hstar.alpha.eta <- -eta * (Hstar + Hstar.eta)

    K <- H / (1 + phi * H)
    Kstar <- Hstar / (1 + phi * Hstar)
    G.phi.phi <- cumsum(c(0, phi * j / (1 + phi * j)^2))[d + 1]

    # derivative of gradient wrt logalpha wrt all parameters except beta
    dl.dlogalpha.dlogalpha <- sum(h.alpha.alpha + Hstar.alpha.alpha /
      (1 + phi * Hstar) - phi * (Hstar.alpha / (1 + phi * Hstar))^2 -
      (1 + phi * d) * (H.alpha.alpha / (1 + phi * H) -
        phi * (H.alpha / (1 + phi * H))^2))
    dl.dlogalpha.dlogeta <- sum(h.alpha.eta + Hstar.alpha.eta /
      (1 + phi * Hstar) - phi * Hstar.alpha * Hstar.eta /
      (1 + phi * Hstar)^2 - (1 + phi * d) * (H.alpha.eta / (1 + phi * H) -
      phi * H.alpha * H.eta / (1 + phi * H)^2))
    dl.dlogalpha.dlogphi <- sum(phi * (-Hstar.alpha * Hstar /
      (1 + phi * Hstar)^2 + H.alpha * H / (1 + phi * H)^2 -
      d * (H.alpha / (1 + phi * H) - phi * H.alpha * H / (1 + phi * H)^2)))
    ddl.dlogalpha <- cbind(
      dl.dlogalpha.dlogalpha, dl.dlogalpha.dlogeta,
      dl.dlogalpha.dlogphi
    )

    # derivative of gradient wrt logeta wrt all parameters except beta
    dl.dlogeta.dlogeta <- sum(h.eta.eta + Hstar.eta.eta / (1 + phi * Hstar) -
      phi * (Hstar.eta / (1 + phi * Hstar))^2 - (1 + phi * d) *
        (H.eta.eta / (1 + phi * H) - phi * (H.eta / (1 + phi * H))^2))
    dl.dlogeta.dlogphi <- sum(phi * (-Hstar.eta * Hstar /
      (1 + phi * Hstar)^2 + H.eta * H / (1 + phi * H)^2 -
      d * (H.eta / (1 + phi * H) - phi * H.eta * H / (1 + phi * H)^2)))
    ddl.dlogeta <- cbind(
      dl.dlogalpha.dlogeta, dl.dlogeta.dlogeta,
      dl.dlogeta.dlogphi
    )

    # derivative of gradient wrt logphi wrt all parameters except beta
    dl.dlogphi.dlogphi <- sum(G.phi.phi + 1 / phi *
      (log(1 + phi * Hstar) - log(1 + phi * H)) + K - Kstar + phi
    * (K^2 - Kstar^2) + d * phi * K * (phi * K - 1))
    ddl.dlogphi <- cbind(
      dl.dlogalpha.dlogphi, dl.dlogeta.dlogphi,
      dl.dlogphi.dlogphi
    )

    # derivatives of gradients wrt (logalpha, logeta, logphi) wrt beta
    H <- (t / alpha)^eta * exp(B)
    Hstar <- (tstar / alpha)^eta * exp(B)
    XX <- c(X) * X[rep(seq_len(nrow(X)), nbeta), ]
    nbeta_rep <- rep(seq_len(nbeta), each = nrow(X))
    H.beta <- as.matrix(aggr(H * X, clusters))
    H.beta2 <- H.beta[rep(seq_len(nrow(H.beta)), nbeta), ] * c(H.beta)
    Hstar.beta <- as.matrix(aggr(Hstar * X, clusters))
    Hstar.beta2 <- Hstar.beta[rep(seq_len(nrow(Hstar.beta)), nbeta), ] * c(Hstar.beta)
    Hstar.beta.beta <- data.table(nbeta_rep, clusters, Hstar * XX)
    H.beta.beta <- data.table(nbeta_rep, clusters, H * XX)
    H <- aggr(H, clusters)
    Hstar <- aggr(Hstar, clusters)
    Hstar.beta.beta <- data.table(clusters, nbeta_rep, Hstar.beta.beta)
    Hstar.beta.beta <- as.matrix(Hstar.beta.beta[,
      j = lapply(.SD, sum),
      by = list(nbeta_rep, clusters)
    ])[, -1:-2, drop = FALSE]
    H.beta.beta <- data.table(clusters, nbeta_rep, H.beta.beta)
    H.beta.beta <- as.matrix(H.beta.beta[,
      j = lapply(.SD, sum),
      by = list(nbeta_rep, clusters)
    ])[, -1:-2, drop = FALSE]
    H.alpha.beta <- -eta * H.beta
    Hstar.alpha.beta <- -eta * Hstar.beta

    dl.dlogalpha.dlogbeta <- colSums(as.matrix(Hstar.alpha.beta /
      (1 + phi * Hstar) - phi * Hstar.alpha * Hstar.beta /
        (1 + phi * Hstar)^2 - (1 + phi * d) * (H.alpha.beta /
        (1 + phi * H) - phi * H.alpha * H.beta / (1 + phi * H)^2)))
    ddl.dlogalpha <- cbind(ddl.dlogalpha, t(dl.dlogalpha.dlogbeta))


    dl.dlogeta.dlogbeta <- t(colSums(as.matrix(Hstar.eta.beta /
      (1 + phi * Hstar) - phi * Hstar.eta * Hstar.beta /
        (1 + phi * Hstar)^2 - (1 + phi * d) * (H.eta.beta / (1 + phi * H) -
        phi * H.eta * H.beta / (1 + phi * H)^2))))
    ddl.dlogeta <- cbind(ddl.dlogeta, dl.dlogeta.dlogbeta)


    dl.dlogphi.dlogbeta <- t(colSums(as.matrix(phi *
      (-Hstar.beta * Hstar / (1 + phi * Hstar)^2 + H.beta * H /
        (1 + phi * H)^2 - d * (H.beta / (1 + phi * H) - phi * H.beta * H /
        (1 + phi * H)^2)))))
    ddl.dlogphi <- cbind(ddl.dlogphi, dl.dlogphi.dlogbeta)

    # derivative of gradient wrt to beta wrt to all parameters
    dl.dlogbeta.dlogbeta <- (Hstar.beta.beta / (1 + phi * Hstar) -
      phi * Hstar.beta2 / (1 + phi * Hstar)^2 - (1 + phi * d) *
        (H.beta.beta / (1 + phi * H) - phi * H.beta2 / (1 + phi * H)^2))
    nbeta_rep2 <- rep(1:nbeta, each = length(H))
    dl.dlogbeta.dlogbeta <- data.table(nbeta_rep2, dl.dlogbeta.dlogbeta)
    dl.dlogbeta.dlogbeta <- as.matrix(dl.dlogbeta.dlogbeta[,
      j = lapply(.SD, sum), by = list(nbeta_rep2)
    ])[, -1]
    ddl.dlogbeta <- rbind(
      dl.dlogalpha.dlogbeta, dl.dlogeta.dlogbeta,
      dl.dlogphi.dlogbeta, dl.dlogbeta.dlogbeta
    )

    hessian <- cbind(
      t(ddl.dlogalpha), t(ddl.dlogeta), t(ddl.dlogphi),
      ddl.dlogbeta
    )
    return(hessian)
  }

  #---PREPARATION---
  call <- match.call()

  # delete rows with missing on variables in the model
  data.temp <- data
  m <- model.matrix(object = formula, data = data.temp)
  data.temp <- data.temp[match(rownames(m), rownames(data.temp)), ]
  X <- model.matrix(formula, data = data.temp)[, -1, drop = FALSE]
  clusters <- data.temp[, clusterid]
  n <- nrow(X)
  ncluster <- length(unique(clusters))
  nbeta <- ncol(X)
  npar <- 3 + nbeta
  if (missing(init)) init <- c(rep(0, npar))

  # extract start variable, end variable and event variable
  Y <- model.extract(
    frame = model.frame(formula = formula, data = data.temp),
    component = "response"
  )
  if (ncol(Y) == 2) {
    tstar <- rep(0, nrow(data.temp))
    t <- Y[, 1]
    delta <- Y[, 2]
  }
  if (ncol(Y) == 3) {
    tstar <- Y[, 1]
    t <- Y[, 2]
    delta <- Y[, 3]
  }

  # number of uncensored in each cluster
  d <- aggr(delta, clusters)
  D <- max(d)
  j <- 0:(D - 1)

  #---MODEL FITTING

  # maximize log likelihood
  fit <- optim(
    par = init, fn = like, gr = gradientfunc, method = "BFGS", hessian = FALSE,
    control = list(fnscale = -1)
  )
  est <- fit$par
  names(est) <- c(
    "log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
    colnames(X)
  )

  # calculate score contributions
  score <- scorefunc(par = est)
  colnames(score) <- c(
    "log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
    colnames(X)
  )

  # calculate hessian
  hessian <- hessianfunc(par = est)
  rownames(hessian) <- c(
    "log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
    colnames(X)
  )
  colnames(hessian) <- c(
    "log(\U003B1)", "log(\U003B7)", "log(\U003D5)",
    colnames(X)
  )

  output <- c(list(
    formula = formula, data = data, clusterid = clusterid,
    ncluster = ncluster, n = n, X = X, fit = fit, est = est,
    score = score, vcov = -solve(hessian), call = call
  ))

  class(output) <- "parfrailty"
  return(output)
}

#' @title Summarizes parfrailty fit
#'
#' @description This is a \code{summary} method for class \code{"parfrailty"}.
#'
#' @param object an object of class \code{"parfrailty"}.
#' @param CI.type string, indicating the type of confidence intervals. Either
#' "plain", which gives untransformed intervals, or "log", which gives
#' log-transformed intervals.
#' @param CI.level desired coverage probability of confidence intervals, in
#' decimal form.
#' @param digits the number of significant digits to use when printing..
#' @param \dots not used.
#' @author Arvid Sjolander and Elisabeth Dahlqwist.
#' @seealso \code{\link{parfrailty}}
#' @examples
#' ## See documentation for parfrailty
#'
#' @rdname summary
#' @export summary.parfrailty
#' @export
summary.parfrailty <- function(object, CI.type = "plain", CI.level = 0.95,
                               digits = max(3L, getOption("digits") - 3L), ...) {
  if (missing(CI.level)) CI.level <- 0.95
  if (missing(CI.type)) CI.type <- "plain"

  ### Inference
  var <- diag(object$vcov)
  se <- sqrt(var)
  zvalue <- object$est / se
  pvalue <- 2 * pnorm(-abs(zvalue))
  confidence.interval <- CI(
    est = object$est, var = var,
    CI.type = CI.type, CI.level = CI.level
  )
  colnames(confidence.interval) <- c("Lower limit", "Upper limit")

  ans <- list(
    est = object$est, se = se, zvalue = zvalue,
    pvalue = pvalue, score = object$score, X = object$X, vcov = object$vcov,
    call = object$call, formula = object$formula, modelmatrix = object$X,
    data = object$data, clusterid = object$clusterid,
    ncluster = object$ncluster, n = object$n, CI.type = CI.type, CI.level = CI.level,
    confidence.interval = confidence.interval
  )

  class(ans) <- "summary.parfrailty"
  return(ans)
}


#' @title Prints summary of parfrailty fit
#'
#' @description This is a \code{print} method for class \code{"summary.parfrailty"}.
#'
#' @param x an object of class \code{"summary.parfrailty"}.
#' @param digits the number of significant digits to use when printing.
#' @param \dots not used.
#' @author Arvid Sjolander and Elisabeth Dahlqwist
#' @seealso \code{\link{parfrailty}}
#' @examples
#'
#' ## See documentation for frailty
#'
#' @rdname print
#' @export print.summary.parfrailty
#' @export
print.summary.parfrailty <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     ...) {
  ## Function call
  cat("Call:", "\n")
  print.default(x$call)
  cat("\nEstimated parameters in the Gamma-Weibull frailty model", "\n")
  cat("\n")
  table.est <- cbind(x$est, exp(x$est), x$se, x$zvalue, x$pvalue)
  rownames(table.est) <- names(x$est)
  colnames(table.est) <- c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")
  # print.default(table.est, digits=3)
  printCoefmat(table.est, digits = 3)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of clusters:", x$ncluster, "\n")
  cat("\n")
}

#' @title Regression standardization in shared frailty gamma-Weibull models
#'
#' @description \code{stdParfrailty} performs regression standardization in shared frailty
#' gamma-Weibull models, at specified values of the exposure, over the sample
#' covariate distribution. Let \eqn{T}, \eqn{X}, and \eqn{Z} be the survival
#' outcome, the exposure, and a vector of covariates, respectively.
#' \code{stdParfrailty} uses a fitted Cox proportional hazards model to
#' estimate the standardized survival function
#' \eqn{\theta(t,x)=E\{S(t|X=x,Z)\}}, where \eqn{t} is a specific value of
#' \eqn{T}, \eqn{x} is a specific value of \eqn{X}, and the expectation is over
#' the marginal distribution of \eqn{Z}.
#'
#' @details \code{stdParfrailty} assumes that a shared frailty gamma-Weibull model
#' \deqn{\lambda(t_{ij}|X_{ij},Z_{ij})=\lambda(t_{ij};\alpha,\eta)U_iexp\{h(X_{ij},Z_{ij};\beta)\}}
#' has been fitted, with parametrization as descibed in the help section for
#' \code{parfrailty}. Integrating out the gamma frailty gives the survival
#' function
#' \deqn{S(t|X,Z)=[1+\phi\Lambda_0(t;\alpha,\eta)exp\{h(X,Z;\beta)\}]^{-1/\phi},}
#' where \eqn{\Lambda_0(t;\alpha,\eta)} is the cumulative baseline hazard
#' \deqn{(t/\alpha)^{\eta}.} The ML estimates of \eqn{(\alpha,\eta,\phi,\beta)}
#' are used to obtain estimates of the survival function \eqn{S(t|X=x,Z)}:
#' \deqn{\hat{S}(t|X=x,Z)=[1+\hat{\phi}\Lambda_0(t;\hat{\alpha},\hat{\eta})exp\{h(X,Z;\hat{\beta})\}]^{-1/\hat{\phi}}.}
#' For each \eqn{t} in the \code{t} argument and for each \eqn{x} in the
#' \code{x} argument, these estimates are averaged across all subjects (i.e.
#' all observed values of \eqn{Z}) to produce estimates
#' \deqn{\hat{\theta}(t,x)=\sum_{i=1}^n \hat{S}(t|X=x,Z_i)/n.} The variance for
#' \eqn{\hat{\theta}(t,x)} is obtained by the sandwich formula.
#'
#' @param fit an object of class \code{"parfrailty"}, as returned by the
#' \code{parfrailty} function in the \pkg{stdReg} package.
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
#' @param clusterid a string containing the name of the cluster identification
#' variable.
#' @param subsetnew an optional logical statement specifying a subset of
#' observations to be used in the standardization. This set is assumed to be a
#' subset of the subset (if any) that was used to fit the regression model.
#' @return An object of class \code{"stdParfrailty"} is a list containing
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
#' \code{stdParfrailty} does not currently handle time-varying exposures or
#' covariates.
#'
#' \code{stdParfrailty} internally loops over all values in the \code{t}
#' argument. Therefore, the function will usually be considerably faster if
#' \code{length(t)} is small.
#'
#' The variance calculation performed by \code{stdParfrailty} does not
#' condition on the observed covariates \eqn{\bar{Z}=(Z_1,...,Z_n)}. To see how
#' this matters, note that
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
#' Dahlqwist E., Pawitan Y., Sjolander A. (2019). Regression standardization
#' and attributable fraction estimation with between-within frailty models for
#' clustered survival data. \emph{Statistical Methods in Medical Research}
#' \bold{28}(2), 462-485.
#'
#' Gail M.H. and Byar D.P. (1986). Variance calculations for direct adjusted
#' survival curves, with applications to testing for no treatement effect.
#' \emph{Biometrical Journal} \bold{28}(5), 587-599.
#'
#' Makuch R.W. (1982). Adjusted survival curve estimation using covariates.
#' \emph{Journal of Chronic Diseases} \bold{35}, 437-443.
#' @examples
#'
#'
#'
#' require(survival)
#'
#' # simulate data
#' set.seed(6)
#' n <- 300
#' m <- 3
#' alpha <- 1.5
#' eta <- 1
#' phi <- 0.5
#' beta <- 1
#' id <- rep(1:n, each = m)
#' U <- rep(rgamma(n, shape = 1 / phi, scale = phi), each = m)
#' X <- rnorm(n * m)
#' # reparametrize scale as in rweibull function
#' weibull.scale <- alpha / (U * exp(beta * X))^(1 / eta)
#' T <- rweibull(n * m, shape = eta, scale = weibull.scale)
#'
#' # right censoring
#' C <- runif(n * m, 0, 10)
#' D <- as.numeric(T < C)
#' T <- pmin(T, C)
#'
#' # strong left-truncation
#' L <- runif(n * m, 0, 2)
#' incl <- T > L
#' incl <- ave(x = incl, id, FUN = sum) == m
#' dd <- data.frame(L, T, D, X, id)
#' dd <- dd[incl, ]
#'
#' fit <- parfrailty(formula = Surv(L, T, D) ~ X, data = dd, clusterid = "id")
#' fit.std <- stdParfrailty(
#'   fit = fit,
#'   data = dd,
#'   X = "X",
#'   x = seq(-1, 1, 0.5),
#'   t = 1:5,
#'   clusterid = "id"
#' )
#' print(summary(fit.std, t = 3))
#' plot(fit.std)
#'
#' @export stdParfrailty
stdParfrailty <- function(fit, data, X, x, t, clusterid, subsetnew) {
  call <- match.call()

  #---PREPARATION---

  formula <- fit$formula
  npar <- length(fit$est)

  # delete rows that did not contribute to the model fit,
  # e.g. missing data or in subset for fit
  # note: object=fit does not work when fit is parfrailty object
  m <- model.matrix(object = formula, data = data)
  data <- data[match(rownames(m), rownames(data)), ]
  n <- nrow(data)

  # Make new subset if supplied.
  subsetnew <-
    if (missing(subsetnew)) {
      rep(1, n)
    } else {
      as.numeric(eval(substitute(subsetnew), data, parent.frame()))
    }
  input <- as.list(environment())

  # extract end variable and event variable
  Y <- model.extract(
    frame = model.frame(formula = formula, data = data),
    component = "response"
  )
  if (ncol(Y) == 2) {
    end <- Y[, 1]
    event <- Y[, 2]
  }
  if (ncol(Y) == 3) {
    end <- Y[, 2]
    event <- Y[, 3]
  }

  clusters <- data[, clusterid]
  ncluster <- length(unique(clusters))

  # assign values to x and reference if not supplied
  # make sure x is a factor if data[, X] is a factor,
  # with the same levels as data[, X]
  if (missing(x)) {
    if (is.factor(data[, X])) {
      x <- as.factor(levels(data[, X]))
    }
    if (is.numeric(data[, X])) {
      if (is.binary(data[, X])) {
        x <- c(0, 1)
      } else {
        x <- round(mean(data[, X], na.rm = TRUE), 2)
      }
    }
  } else {
    if (is.factor(x)) {
      temp <- x
      levels(x) <- levels(data[, X])
      x[seq_len(length(x))] <- temp
    } else {
      if (is.factor(data[, X])) {
        x <- factor(x)
        temp <- x
        levels(x) <- levels(data[, X])
        x[seq_len(length(x))] <- temp
      }
    }
  }
  input$x <- x
  nX <- length(x)

  # assign value to t if missing
  if (missing(t)) {
    stop("You have to specify the times (t) at which to estimate the standardized survival function")
    # t <- end[event == 1]
  }
  t <- sort(t)
  input$t <- t
  nt <- length(t)

  # preparation
  est <- matrix(nrow = nt, ncol = nX)
  vcov <- vector(mode = "list", length = nt)
  logalpha <- fit$est[1]
  alpha <- exp(logalpha)
  logeta <- fit$est[2]
  eta <- exp(logeta)
  logphi <- fit$est[3]
  phi <- exp(logphi)
  beta <- fit$est[(3 + 1):npar]

  #---LOOP OVER nt

  for (j in 1:nt) {
    if (t[j] == 0) {
      est[j, ] <- 1
      vcov[[j]] <- matrix(0, nrow = nX, ncol = nX)
    } else {
      H0t <- (t[j] / alpha)^eta

      #---ESTIMATES OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---

      si <- matrix(nrow = n, ncol = nX)
      SI.logalpha <- vector(length = nX)
      SI.logeta <- vector(length = nX)
      SI.logphi <- vector(length = nX)
      SI.beta <- matrix(nrow = nX, ncol = npar - 3)
      for (i in 1:nX) {
        data.x <- data
        if (!is.na(x[i])) {
          data.x[, X] <- x[i]
        }
        m <- model.matrix(object = formula, data = data.x)[, -1, drop = FALSE]
        predX <- colSums(beta * t(m))
        temp <- 1 + phi * H0t * exp(predX)
        si[, i] <- temp^(-1 / phi)
        SI.logalpha[i] <- mean(subsetnew * H0t * eta * exp(predX) / temp^(1 / phi + 1)) * n /
          ncluster
        SI.logeta[i] <- mean(subsetnew * (-H0t) * exp(predX) * log(t[j] / alpha) * eta /
          temp^(1 / phi + 1)) * n / ncluster
        SI.logphi[i] <- mean(subsetnew * log(temp) / (phi * temp^(1 / phi)) -
          H0t * exp(predX) / temp^(1 / phi + 1)) * n / ncluster
        SI.beta[i, ] <- colMeans(subsetnew * (-H0t) * exp(predX) * m / temp^(1 / phi + 1)) *
          n / ncluster
      }
      est[j, ] <- colSums(subsetnew * si, na.rm = TRUE) / sum(subsetnew)

      #---VARIANCE OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---

      sres <- subsetnew * (si - matrix(rep(est[j, ], each = n), nrow = n, ncol = nX))
      sres <- aggr(sres, clusters)
      coefres <- fit$score
      res <- cbind(sres, coefres)

      J <- var(res, na.rm = TRUE)

      # Note: the term n/ncluster is because SI.logalpha, SI.logeta, SI.logphi,
      # and SI.beta are clustered, which they are not in stdCoxph
      SI <- cbind(
        -diag(nX) * mean(subsetnew) * n / ncluster, SI.logalpha, SI.logeta,
        SI.logphi, SI.beta
      )

      betaI <- cbind(matrix(0, nrow = npar, ncol = nX), -solve(fit$vcov) / ncluster)
      I <- rbind(SI, betaI)
      V <- (solve(I) %*% J %*% t(solve(I)) / ncluster)[1:nX, 1:nX]
      vcov[[j]] <- V
    }
  }
  out <- list(call = call, input = input, est = est, vcov = vcov)
  #---OUTPUT---

  class(out) <- "stdParfrailty"
  return(out)
}

#' @title Summarizes Frailty standardization fit
#'
#' @description This is a \code{summary} method for class \code{"stdParfrailty"}.
#'
#' @param object an object of class \code{"stdParfrailty"}.
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
#' @seealso \code{\link{stdParfrailty}}
#' @examples
#'
#' ## See documentation for stdParfrailty
#'
#' @rdname summary
#' @export summary.stdParfrailty
#' @export
summary.stdParfrailty <- summary.stdCoxph

#' @rdname print
#' @export print.stdParfrailty
#' @export
print.stdParfrailty <- function(x, ...) {
  print(summary(x))
}

#' @title Prints summary of Frailty standardization fit
#'
#' @description This is a \code{print} method for class \code{"summary.stdParfrailty"}.
#'
#'
#' @param x an object of class \code{"summary.stdParfrailty"}.
#' @param \dots not used.
#' @author Arvid Sjolander
#' @seealso \code{\link{stdParfrailty}}
#' @examples
#'
#'
#' ## See documentation for stdParfrailty
#'
#' @rdname print
#' @export print.summary.stdParfrailty
#' @export
print.summary.stdParfrailty <- print.summary.stdCoxph

#' @title Plots parfrailty standardization fit
#'
#' @description This is a \code{plot} method for class \code{"stdParfrailty"}.
#'
#' @param x an object of class \code{"stdParfrailty"}.
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
#' @seealso \code{\link{stdParfrailty}}
#' @examples
#'
#'
#' ## See documentation for stdParfrailty
#'
#' @rdname plot
#' @export plot.stdParfrailty
#' @export
plot.stdParfrailty <- plot.stdCoxph
