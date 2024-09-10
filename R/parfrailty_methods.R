#' @title Fits shared frailty gamma-Weibull models
#'
#' @description \code{parfrailty} fits shared frailty gamma-Weibull models. It is
#' specifically designed to work with the function \code{standardize_parfrailty}, which
#' performs regression standardization in shared frailty gamma-Weibull models.
#'
#' @details \code{parfrailty} fits the shared frailty gamma-Weibull model
#' \deqn{\lambda(t_{ij}|C_{ij})=\lambda(t_{ij};\alpha,\eta)U_i\exp\{h(C_{ij};\beta)\},}
#' where \eqn{t_{ij}} and \eqn{C_{ij}} are the survival time and covariate
#' vector for subject \eqn{j} in cluster \eqn{i}, respectively.
#' \eqn{\lambda(t;\alpha,\eta)} is the Weibull baseline hazard function
#' \deqn{\eta t^{\eta-1}\alpha^{-\eta},} where \eqn{\eta} is the shape
#' parameter and \eqn{\alpha} is the scale parameter. \eqn{U_i} is the
#' unobserved frailty term for cluster \eqn{i}, which is assumed to have a
#' gamma distribution with scale = 1/shape = \eqn{\phi}. \eqn{h(X;\beta)} is
#' the regression function as specified by the \code{formula} argument,
#' parameterized by a vector \eqn{\beta}. The ML estimates
#' \eqn{\{\log(\hat{\alpha}),\log(\hat{\eta}),\log(\hat{\phi}),\hat{\beta}\}} are
#' obtained by maximizing the marginal (over \eqn{U}) likelihood.
#'
#' @param formula an object of class "\code{formula}", in the same format as
#' accepted by the \link[survival]{coxph} function.
#' @param data a data frame containing the variables in the model.
#' @param clusterid a string containing the name of a cluster identification
#' variable.
#' @param init an optional vector of initial values for the model parameters.
#' @returns An object of class \code{"parfrailty"} which is a list containing:
#' \item{est}{ the Maximum Likelihood (ML) estimates \eqn{\{\log(\hat{\alpha}),\log(\hat{\eta}),
#' \log(\hat{\phi}),\hat{\beta}\}}. } \item{vcov}{ the variance-covariance
#' vector of the ML estimates. } \item{score}{ a matrix containing the
#' cluster-specific contributions to the ML score equations.  }
#' @note If left truncation is present, it is assumed that it is strong left
#' truncation.  This means that even if the truncation time may be
#' subject-specific, the whole cluster is unobserved if at least one subject in
#' the cluster dies before his/her truncation time. If all subjects in the
#' cluster survive beyond their subject-specific truncation times, then the
#' whole cluster is observed (Van den Berg and Drepper, 2016).
#' @author Arvid Sjölander and Elisabeth Dahlqwist.
#' @references Dahlqwist E., Pawitan Y., Sjölander A. (2019). Regression
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
#' # reparameterize scale as in rweibull function
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
#' print(fit)
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
    h <- delta * log(eta * times^(eta - 1) / alpha^eta * exp(B))
    H <- (times / alpha)^eta * exp(B)
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
    h.eta <- delta * (1 + eta * (log(times) - log(alpha)))
    h.beta <- X * delta
    H <- (times / alpha)^eta * exp(B)
    Hstar <- (tstar / alpha)^eta * exp(B)
    H.eta <- eta * log(times / alpha) * H
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
    XX <- c(X) * X[rep(seq_len(nrow(X)), nbeta), ]
    h.eta <- delta * (1 + eta * (log(times) - log(alpha)))
    H <- (times / alpha)^eta * exp(B)
    Hstar <- (tstar / alpha)^eta * exp(B)
    H.eta <- eta * log(times / alpha) * H
    Hstar.eta <- eta * log(tstar / alpha) * Hstar
    Hstar.eta[tstar == 0] <- 0
    H.eta.eta <- H.eta + eta^2 * (log(times / alpha))^2 * H
    Hstar.eta.eta <- Hstar.eta + eta^2 * (log(tstar / alpha))^2 * Hstar
    Hstar.eta.eta[tstar == 0] <- 0
    H.eta.beta <- eta * log(times / alpha) * (H * X)
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
    H <- (times / alpha)^eta * exp(B)
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
    times <- Y[, 1]
    delta <- Y[, 2]
  }
  if (ncol(Y) == 3) {
    tstar <- Y[, 1]
    times <- Y[, 2]
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
#' @param ci_type string, indicating the type of confidence intervals. Either
#' "plain", which gives untransformed intervals, or "log", which gives
#' log-transformed intervals.
#' @param ci_level desired coverage probability of confidence intervals, in
#' decimal form.
#' @param digits the number of significant digits to use when printing..
#' @param \dots not used.
#' @author Arvid Sjölander and Elisabeth Dahlqwist.
#' @seealso \code{\link{parfrailty}}
#' @examples
#' ## See documentation for parfrailty
#'
#' @rdname summary
#' @export summary.parfrailty
#' @export
#' @returns An object of class "summary.parfrailty", which is a list that contains relevant summary statistics about the fitted model
summary.parfrailty <- function(object, ci_type = "plain", ci_level = 0.95,
                               digits = max(3L, getOption("digits") - 3L), ...) {
  if (missing(ci_level)) ci_level <- 0.95
  if (missing(ci_type)) ci_type <- "plain"

  ### Inference
  var <- diag(object$vcov)
  se <- sqrt(var)
  zvalue <- object$est / se
  pvalue <- 2 * pnorm(-abs(zvalue))
  confidence.interval <- CI(
    est = object$est, var = var,
    ci_type = ci_type, ci_level = ci_level
  )
  colnames(confidence.interval) <- c("Lower limit", "Upper limit")

  ans <- list(
    est = object$est, se = se, zvalue = zvalue,
    pvalue = pvalue, score = object$score, X = object$X, vcov = object$vcov,
    call = object$call, formula = object$formula, modelmatrix = object$X,
    data = object$data, clusterid = object$clusterid,
    ncluster = object$ncluster, n = object$n, ci_type = ci_type, ci_level = ci_level,
    confidence.interval = confidence.interval
  )

  class(ans) <- "summary.parfrailty"
  return(ans)
}

#' Print method for parametric frailty fits
#'
#' @param x An object of class "parfrailty"
#' @param digits Number of digits to print
#' @param ... Not used
#' @export
#' @returns The object being printed, invisibly
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
  printCoefmat(table.est, digits = 3)
  cat("\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of clusters:", x$ncluster, "\n")
  cat("\n")
  invisible(x)
}

#' @title Regression standardization in shared frailty gamma-Weibull models
#'
#' @description \code{standardize_parfrailty} performs regression standardization in shared frailty
#' gamma-Weibull models, at specified values of the exposure, over the sample
#' covariate distribution. Let \eqn{T}, \eqn{X}, and \eqn{Z} be the survival
#' outcome, the exposure, and a vector of covariates, respectively.
#' \code{standardize_parfrailty} fits a parametric frailty model to
#' estimate the standardized survival function
#' \eqn{\theta(t,x)=E\{S(t|X=x,Z)\}}, where \eqn{t} is a specific value of
#' \eqn{T}, \eqn{x} is a specific value of \eqn{X}, and the expectation is over
#' the marginal distribution of \eqn{Z}.
#'
#' @details \code{standardize_parfrailty} fits a shared frailty gamma-Weibull model
#' \deqn{\lambda(t_{ij}|X_{ij},Z_{ij})=\lambda(t_{ij};\alpha,\eta)U_iexp\{h(X_{ij},Z_{ij};\beta)\}}
#' , with parameterization as described in the help section for
#' \link{parfrailty}. Integrating out the gamma frailty gives the survival
#' function
#' \deqn{S(t|X,Z)=[1+\phi\Lambda_0(t;\alpha,\eta)\exp\{h(X,Z;\beta)\}]^{-1/\phi},}
#' where \eqn{\Lambda_0(t;\alpha,\eta)} is the cumulative baseline hazard
#' \deqn{(t/\alpha)^{\eta}.} The ML estimates of \eqn{(\alpha,\eta,\phi,\beta)}
#' are used to obtain estimates of the survival function \eqn{S(t|X=x,Z)}:
#' \deqn{\hat{S}(t|X=x,Z)=[1+\hat{\phi}\Lambda_0(t;\hat{\alpha},\hat{\eta})\exp\{h(X,Z;\hat{\beta})\}]^{-1/\hat{\phi}}.}
#' For each \eqn{t} in the \code{t} argument and for each \eqn{x} in the
#' \code{x} argument, these estimates are averaged across all subjects (i.e.
#' all observed values of \eqn{Z}) to produce estimates
#' \deqn{\hat{\theta}(t,x)=\sum_{i=1}^n \hat{S}(t|X=x,Z_i)/n.} The variance for
#' \eqn{\hat{\theta}(t,x)} is obtained by the sandwich formula.
#'
#' @inherit standardize_coxph
#' @author Arvid Sjölander
#' @references
#'
#' Chang I.M., Gelman G., Pagano M. (1982). Corrected group prognostic curves
#' and summary statistics. \emph{Journal of Chronic Diseases} \bold{35},
#' 669-674.
#'
#' Dahlqwist E., Pawitan Y., Sjölander A. (2019). Regression standardization
#' and attributable fraction estimation with between-within frailty models for
#' clustered survival data. \emph{Statistical Methods in Medical Research}
#' \bold{28}(2), 462-485.
#'
#' Gail M.H. and Byar D.P. (1986). Variance calculations for direct adjusted
#' survival curves, with applications to testing for no treatment effect.
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
#' # reparameterize scale as in rweibull function
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
#' fit.std <- standardize_parfrailty(
#'   formula = Surv(L, T, D) ~ X,
#'   data = dd,
#'   values = list(X = seq(-1, 1, 0.5)),
#'   times = 1:5,
#'   clusterid = "id"
#' )
#' print(fit.std)
#' plot(fit.std)
#'
#' @export standardize_parfrailty
standardize_parfrailty <- function(formula,
                                   data,
                                   values,
                                   times,
                                   clusterid,
                                   ci_level = 0.95,
                                   ci_type = "plain",
                                   contrasts = NULL,
                                   family = "gaussian",
                                   reference = NULL,
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
      parfrailty(formula = formula, data = data, clusterid = clusterid)
    },
    error = function(cond) {
      return(cond)
    }
  )
  if (inherits(fit, "simpleError")) {
    stop("parfrailty function failed with error: ", fit[["message"]])
  }

  #---PREPARATION---
  npar <- length(fit$est)

  # delete rows that did not contribute to the model fit,
  # e.g. missing data or in subset for fit
  # note: object=fit does not work when fit is parfrailty object
  m <- model.matrix(object = formula, data = data)
  data <- data[match(rownames(m), rownames(data)), ]
  n <- nrow(data)

  input <- as.list(environment())

  clusters <- data[, clusterid]
  ncluster <- length(unique(clusters))

  nX <- nrow(valuesout)

  # assign value to times if missing
  if (missing(times)) {
    stop("You have to specify the times at which to estimate the standardized survival function")
    # times <- end[event == 1]
  }
  times <- sort(times)
  input$times <- times
  nt <- length(times)

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
  si <- array(dim = c(n,nX,nt))
  SI.logalpha <- matrix(nrow = nX, ncol = nt)
  SI.logeta <- matrix(nrow = nX, ncol = nt)
  SI.logphi <- matrix(nrow = nX, ncol = nt)
  SI.beta <- array(dim = c(nX, npar - 3, nt))

  for (j in 1:nt) {
    if (times[j] == 0) {
      est[j, ] <- 1
      vcov[[j]] <- matrix(0, nrow = nX, ncol = nX)
    } else {
      H0t <- (times[j] / alpha)^eta

      #---ESTIMATES OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---

      for (i in 1:nX) {
        data.x <- do.call("transform", c(
          list(data),
          valuesout[i, , drop = FALSE]
        ))
        m <- model.matrix(object = formula, data = data.x)[, -1, drop = FALSE]
        predX <- colSums(beta * t(m))
        temp <- 1 + phi * H0t * exp(predX)
        si[, i, j] <- temp^(-1 / phi)
        SI.logalpha[i, j] <- mean(H0t * eta * exp(predX) / temp^(1 / phi + 1)) * n /
          ncluster
        SI.logeta[i, j] <- mean((-H0t) * exp(predX) * log(times[j] / alpha) * eta /
          temp^(1 / phi + 1)) * n / ncluster
        SI.logphi[i, j] <- mean(log(temp) / (phi * temp^(1 / phi)) -
          H0t * exp(predX) / temp^(1 / phi + 1)) * n / ncluster
        SI.beta[i, ,j] <- colMeans((-H0t) * exp(predX) * m / temp^(1 / phi + 1)) *
          n / ncluster
      }
      est[j, ] <- colSums(si[,,j], na.rm = TRUE) / n

      #---VARIANCE OF SURVIVAL PROBABILITIES AT VALUES SPECIFIED BY x ---
    }
  }
  sres <- data.table::rbindlist(lapply(1:nt, function(j) data.table(val = si[,,j] - matrix(rep(est[j, ], each = n),
                                                                               nrow = n, ncol = nX),
                                                        times = j, clusters = clusters)))
  sres <- sres[, j = lapply(.SD, sum), by = list(clusters,times)][,c(-1,-2)]

  n_unique_clusters <- length(unique(clusters))
  for (j in 1:nt){
    if (times[j] !=0){
      coefres <- fit$score
      res <- cbind(sres[((j-1)*n_unique_clusters+1):(j*n_unique_clusters),], coefres)

      J <- var(res, na.rm = TRUE)

      # Note: the term n/ncluster is because SI.logalpha, SI.logeta, SI.logphi,
      # and SI.beta are clustered, which they are not in stdCoxph
      SI <- cbind(
        -diag(nX) * n / ncluster, SI.logalpha[, j], SI.logeta[, j],
        SI.logphi[, j],SI.beta[,, j]
      )

      betaI <- cbind(matrix(0, nrow = npar, ncol = nX), -solve(fit$vcov) / ncluster)
      I <- rbind(SI, betaI)
      V <- (solve(I) %*% J %*% t(solve(I)) / ncluster)[1:nX, 1:nX]
      vcov[[j]] <- V
    }
  }
  out <- list(call = call, input = input, measure = "survival",
              est = est, vcov = vcov)
  #---OUTPUT---

  class(out) <- "std_surv"
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
