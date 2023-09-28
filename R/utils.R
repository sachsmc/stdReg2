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

## copy paste from sandwich package; we don't need the entire package
estfun_glm <- function(x, ...)
{
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) xmat <- xmat[, !alias, drop = FALSE]
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  dispersion <- if(substr(x$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
  else sum(wres^2, na.rm = TRUE)/sum(weights(x, "working"), na.rm = TRUE)
  rval <- wres * xmat / dispersion
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  res <- residuals(x, type = "pearson")
  return(rval)
}

aggr <- function(x, clusters) {
  temp <- data.table(x)
  temp <- as.matrix(temp[, j = lapply(.SD, sum), by = clusters])[, -1]
}

CI <- function(est, var, ci_type = "plain", ci_level = 0.95) {
  se <- sqrt(var)
  qqq <- abs(qnorm((1 - ci_level) / 2))

  if (ci_type == "plain") {
    lower <- est - qqq * se
    upper <- est + qqq * se
  }
  if (ci_type == "log") {
    lower <- est * exp(-qqq * se / est)
    upper <- est * exp(qqq * se / est)
  }

  ci <- cbind(lower, upper)
  return(ci)
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

check_values_data <- function(values, data) {
  xnms <- names(values)
  fitnms <- names(data)
  mnams <- match(xnms, fitnms)
  if (anyNA(mnams)) {
    stop(
      "variable(s) ", toString(xnms[which(is.na(mnams))]),
      " not found in ", deparse1(substitute(fit_outcome)), "$data"
    )
  }
}

get_outcome_exposure <- function(formula_outcome, data, values) {
  outcome <- data[, as.character(formula_outcome)[[2L]]]
  if (!is.data.frame(values)) {
    valuesout <- expand.grid(values)
  } else {
    valuesout <- values
  }
  exposure_names <- colnames(valuesout)
  exposure <- data[, exposure_names]
  list(outcome = outcome, exposure = exposure, exposure_names = exposure_names, valuesout = valuesout)
}

format_result_standardize <- function(res,
                                      contrasts,
                                      references,
                                      transforms,
                                      ci_type,
                                      ci_level,
                                      format_class,
                                      summary_fun_name) {
  contrast <- reference <- NULL
  ## change contrasts, references and transforms to NULL in string format
  if (is.null(contrasts) && !is.null(references) || !is.null(contrasts) && is.null(references)) {
    warning("Reference level or contrast not specified. Defaulting to NULL. ")
  }

  contrasts <- unique(c("NULL", contrasts))
  references <- unique(c("NULL", references))
  transforms <- unique(c("NULL", transforms))
  grid <- expand.grid(
    contrast = contrasts,
    reference = references,
    transform = transforms
  )
  grid <- subset(grid, (contrast == "NULL" & reference == "NULL") |
                   (contrast != "NULL" & reference != "NULL"))
  summary_fun <- function(contrast, reference, transform) {
    null_helper <- function(x) {
      if (is.null(x) || x == "NULL") {
        NULL
      } else {
        as.character(x)
      }
    }
    transform <- null_helper(transform)
    contrast <- null_helper(contrast)
    reference <- null_helper(reference)

    do.call(summary_fun_name, list(object=res,
                                   ci_type = ci_type,
                                   ci_level = ci_level,
                                   transform = transform,
                                   contrast = contrast,
                                   reference = reference
    ))
  }
  res_contrast <- as.list(as.data.frame(do.call(mapply, c("summary_fun", unname(as.list(grid))))))
  res_fin <- list(res_contrast = res_contrast, res = res)
  class(res_fin) <- format_class
  res_fin
}
