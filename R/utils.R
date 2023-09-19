#' Compute the sandwich variance components from a model fit
#'
#' @param fit A fitted model object of class glm, coxph, ah, or survfit
#' @param data The data used to fit the model
#' @param weights Optional weights
#' @param t Optional fixed time point for survival objects
#' @param fit.detail Additional information for survival objects, see Details
#'
#' @return A list consisting of the Fisher information matrix (I) and the Score equations (U)
#' @export sandwich
sandwich <- function(fit, data, weights, t, fit.detail) {
  n <- nrow(data)
  if (missing(weights)) {
    weights <- rep(1, n)
  }

  if (inherits(x = fit, what = "glm")) {
    ## Dispersion parameter = 1. Why?
    m <- expand(model.matrix(fit), rownames(data))
    res <- expand(residuals(fit, type = "response"), rownames(data))
    U <- weights * m * res
    U[is.na(U)] <- 0
    ## Derive Fisher information matrix from asymptotic covariance matrix (MLE theory)
    ## NOTE: summary(fit)$cov.unscaled is weighted
    I <- -solve(summary(fit)$cov.unscaled) / n
  }
  if (inherits(x = fit, what = "ah")) {
    #---meat---

    # residuals are weighted
    res <- predict(object = fit, type = "residuals")
    rownames(res) <- fit$incl
    colnames(res) <- names(fit$coefficients)
    res <- expand(res, rownames(data))
    U <- res
  }
  if (inherits(x = fit, what = "coxph")) {
    #---meat for regression coefficients---

    # score residuals are unweighted, but computed with estimated coefficients
    # from weighted model
    res <- residuals(fit, type = "score")
    res <- expand(res, rownames(data))
    Ucoef <- weights * res

    #---bread for regression coefficients---

    # vcov(fit) is weighted
    Icoef <- -solve(vcov(fit)) / n

    if (missing(t)) {
      U <- Ucoef
      colnames(U) <- names(fit$coef)
      U[is.na(U)] <- 0
      I <- Icoef
    } else {
      nt <- length(t)

      #---meat and bread for baseline hazard---

      # check if left truncation
      varsLHS <- all.vars(fit$formula[[2]])
      t2 <- varsLHS[length(varsLHS) - 1]
      if (length(varsLHS) == 3) {
        t1 <- varsLHS[1]
      } else {
        t1 <- NULL
      }
      time <- fit.detail$time
      # nevent is unweighted
      nevent <- fit.detail$nevent
      # hazard is weighted
      dH <- fit.detail$hazard
      names(dH) <- time
      # varhaz is dH/mean(exp(x*b)) where mean is in risk set
      # varhaz is weighted
      dHvar <- fit.detail$varhaz
      Hvar <- stepfun(time, c(0, cumsum(dHvar)))
      p <- predict(object = fit, type = "risk")
      m <- model.matrix(fit)
      names(p) <- rownames(m)
      p <- expand(p, rownames(data))
      ncoef <- length(fit$coef)
      # fit.detail$means is weighted, but not relative to the mean covariate
      # in the sample, like all other outputs from coxph.detail,
      # subtracting fit$means fixes this
      means <- as.matrix(fit.detail$means)
      means <- means - matrix(fit$means,
        nrow = nrow(means), ncol = ncol(means),
        byrow = TRUE
      )
      means <- expand(means, data[, t2])
      UH <- matrix(nrow = n, ncol = nt)
      IH <- matrix(nrow = nt, ncol = ncoef)
      for (j in 1:nt) {
        # dividing with nevent accounts for ties,
        # but this assumes that H is computed with Breslow method for ties,
        # and that weights are equal within ties.
        tmp1 <- n * expand(dH / nevent * (time <= t[j]), data[, t2])
        tmp1[is.na(tmp1)] <- 0
        tmp2 <- n * Hvar(pmin(t[j], data[, t2])) * p
        if (!is.null(t1)) {
          tmp2 <- tmp2 - n * Hvar(data[, t1]) * (data[, t1] < t[j]) * p
        }
        UH[, j] <- tmp1 - tmp2
        dH.dbeta <- means * tmp1
        dH.dbeta[is.na(dH.dbeta)] <- 0
        IH[j, ] <- -colMeans(dH.dbeta)
      }
      U <- cbind(Ucoef, UH)
      colnames(U) <- c(names(fit$coef), paste0("H", t))
      U[is.na(U)] <- 0
      I <- rbind(
        cbind(Icoef, matrix(0, nrow = ncoef, ncol = length(t))),
        cbind(IH, -diag(length(t)))
      )
    }
  }
  if (inherits(x = fit, what = "survfit")) {
    #---meat---

    # check if left truncation
    # survfit object has no formula element, so need to get it from call,
    # need to use eval, since the fit$call$formula will be literary what the user
    # gave as argument, e.g. if formula=f, then fit$call$formula is f, not the
    # formula contained in f
    varsLHS <- all.vars(eval(fit$call$formula)[[2]])
    t2 <- varsLHS[length(varsLHS) - 1]
    if (length(varsLHS) == 3) {
      t1 <- varsLHS[1]
    } else {
      t1 <- NULL
    }
    # need to use summary(fit), since n.events and n.risk from fit
    # behave strange when there is left truncation
    ss <- summary(fit)
    strata <- ss$strata
    # n.strata is unweighted
    n.strata <- summary(strata)
    K <- length(n.strata)
    names.strata <- names(n.strata)
    time <- ss$time
    # n.event and n.risk are weighted
    n.event <- ss$n.event
    n.risk <- ss$n.risk
    dH <- n.event / n.risk
    names(dH) <- paste(time, strata)
    dHvar <- dH / n.risk
    # survfit object has no formula element, so need to get it from call,
    # need to use eval, since the fit$call$formula will be literary what the user
    # gave as argument, e.g. if formula=f, then fit$call$formula is f, not the
    # formula contained in f
    vars <- all.vars(eval(fit$call$formula)[[3]])
    # note: strata is a function in the survival package
    strata.all <- strata(data[, vars, drop = FALSE])
    tmp1 <- matrix(nrow = n, ncol = K)
    U <- matrix(nrow = n, ncol = K)
    colnames(U) <- names.strata
    breaks <- c(0, cumsum(n.strata))
    for (k in 1:K) {
      incl <- (breaks[k] + 1):breaks[k + 1]
      Hvar <- stepfun(time[incl], c(0, cumsum(dHvar[incl])))
      # dividing with nevent[incl] account for ties,
      # but this assumes that H is computed with Breslow method for ties,
      # and that weights are equal within ties.
      # multiplying with weights corrects for nevent being weighted;
      # here we just want to divide with the actual number of events to account
      # for ties, not the weighted number of events
      tmp1.time <- n * dH[incl] / n.event[incl] * (time[incl] <= t)
      tmp1[, k] <- tmp1.time[match(
        paste(data[, t2], strata.all),
        names(tmp1.time)
      )] * weights
      tmp1[is.na(tmp1[, k]), k] <- 0
      sk <- names.strata[k]
      incl <- which(strata.all == sk)
      tmp2 <- n * Hvar(pmin(t, data[incl, t2]))
      if (!is.null(t1)) {
        tmp2 <- tmp2 - n * Hvar(data[incl, t1]) * (data[incl, t1] < t)
      }
      U[incl, k] <- tmp1[incl, k] - tmp2
    }

    #---bread---

    I <- diag(-1, K)
    rownames(I) <- names.strata
    colnames(I) <- names.strata
  }

  U[is.na(U)] <- 0
  return(list(I = I, U = U))
}

aggr <- function(x, clusters) {
  temp <- data.table(x)
  temp <- as.matrix(temp[, j = lapply(.SD, sum), by = clusters])[, -1]
}

confintall <- function(object, parm, level = 0.95, fun, type = "plain", ...) {
  est <- do.call(what = fun, args = list(est = object$est))
  var <- delmet(fun = fun, est = object$est, vcov = object$vcov)
  ci <- CI(est = est, var = var, CI.type = type, CI.level = level)
  return(ci)
}

CI <- function(est, var, CI.type = "plain", CI.level = 0.95) {
  se <- sqrt(var)
  qqq <- abs(qnorm((1 - CI.level) / 2))

  if (CI.type == "plain") {
    lower <- est - qqq * se
    upper <- est + qqq * se
  }
  if (CI.type == "log") {
    lower <- est * exp(-qqq * se / est)
    upper <- est * exp(qqq * se / est)
  }

  ci <- cbind(lower, upper)
  return(ci)
}

confint.stdCoxph <- confintall
confint.stdGee <- confintall
confint.stdGlm <- confintall
confint.stdParfrailty <- confintall

delmet <- function(fun, est, vcov) {
  if (!is.list(vcov)) {
    est <- matrix(est, nrow = 1)
    vcov <- list(vcov)
  }
  p <- length(vcov)
  var <- vector(length = p)
  for (i in 1:p) {
    gradient <- matrix(grad(func = fun, x = est[i, , drop = FALSE]), nrow = 1)
    var[i] <- gradient %*% vcov[[i]] %*% t(gradient)
  }
  return(var)
}

expand <- function(x, names) {
  if (is.vector(x)) {
    x <- x[match(names, names(x))]
  }
  if (is.matrix(x)) {
    x <- x[match(names, rownames(x)), , drop = FALSE]
  }
  return(x)
}

Hfun <- function(fit, data, fit.detail) {
  if (inherits(x = fit, what = "survfit")) {
    # need to use summary(fit), since n.events and n.risk from fit
    # behave strange when there is left truncation
    ss <- summary(fit)
    strata <- ss$strata
    # n.strata is unweighted
    n.strata <- summary(strata)
    K <- length(n.strata)
    names.strata <- names(n.strata)
    time <- ss$time
    # n.event and n.risk are weighted
    n.event <- ss$n.event
    n.risk <- ss$n.risk
    dH <- n.event / n.risk
    H <- list()
    breaks <- c(0, cumsum(n.strata))
    for (k in 1:K) {
      incl <- (breaks[k] + 1):breaks[k + 1]
      H[[k]] <- stepfun(time[incl], c(0, cumsum(dH[incl])))
    }
    names(H) <- names.strata
  }
  if (inherits(x = fit, what = "coxph")) {
    time <- fit.detail$time
    # dH is weighted
    dH <- fit.detail$hazard
    H <- stepfun(time, c(0, cumsum(dH)))
  }

  return(H)
}

is.binary <- function(v) {
  if (is.numeric(v) && all(v == 0 | v == 1, na.rm = TRUE)) {
    TRUE
  } else {
    FALSE
  }
}

logit <- function(x) log(x) - log(1 - x)

odds <- function(x) x / (1 - x)
