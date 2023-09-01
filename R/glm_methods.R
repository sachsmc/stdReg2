#' Get regression standardized estimates from a glm
#'
#' @param fit The fitted object as returned by glm
#' @param values A named list or data.frame specifying the variables and values
#'                at which marginal means of the outcome will be estimated.
#' @param weights Optional vector of weights used in the standardization process
#'
#' @returns An object of class std. This is basically a list with components estimates and covariance
#' @examples
#'
#' n <- 1000
#' Z <- rnorm(n)
#' X <- rnorm(n, mean=Z)
#' Y <- rbinom(n, 1, prob=(1+exp(X+Z))^(-1))
#' dd <- data.frame(Z, X, Y)
#' fit <- glm(formula=Y~X+Z+X*Z, family="binomial", data=dd)
#' standardize.glm(fit=fit, values = list(X = 0:1))
#'
standardize.glm <- function(fit, values, weights) {

  if(is.null(fit$data)) {
    stop("data element missing from fit")
  }
  xnms <- names(values)
  fitnms <- names(fit$data)
  mnams <- match(xnms, fitnms)
  if(any(is.na(mnams))) {
    stop(paste0("variable(s) ", paste(xnms[which(is.na(mnams))], collapse = ", "),
                " not found in ", deparse1(substitute(fit)), "$data"))
  }

  n <- nrow(fit$data)

  if(missing(weights)) {
    weights <- rep(1, n)
  }

  if(!is.data.frame(values)) {
    valuesout <- expand.grid(values)
  } else valuesout <- values

  predmat <- matrix(nrow = n, ncol = nrow(valuesout))
  g <- family(fit)$mu.eta
  SI.beta <- matrix(nrow=nrow(valuesout), ncol=length(fit$coefficients))

  for(i in 1:nrow(valuesout)){
    data.x <- do.call("transform", c(list(fit$data),
                                     valuesout[i,,drop = FALSE]))

    pred.x <- predict(object = fit, newdata = data.x)
    predmat[, i] <- fit$family$linkinv(pred.x)

    dmu.deta <- g(pred.x)

    deta.dbeta <- model.matrix(object=terms(fit), data=data.x)
    dmu.dbeta <- dmu.deta*deta.dbeta
    SI.beta[i, ] <- colMeans(weights*dmu.dbeta)
  }
  est <- colSums(weights*predmat, na.rm=TRUE)/sum(weights)

  sandwich.fit <- sandwich(fit=fit, data=fit$data, weights=weights)
  ores <- sandwich.fit$U
  mres <- weights * (predmat - matrix(rep(est, each=n), nrow=n,
                                         ncol=nrow(valuesout)))
  res <- cbind(mres, ores)

  SI <- cbind(-diag(nrow(valuesout))*mean(weights), SI.beta)
  oI <- cbind(matrix(0, nrow=length(fit$coefficients),
                     ncol=nrow(valuesout)), sandwich.fit$I)
  J <- var(res, na.rm=TRUE)
  I <- rbind(SI, oI)
  variance <- (solve(I)%*%J%*%t(solve(I))/n)[1:nrow(valuesout), 1:nrow(valuesout)]

  rownames(variance) <- colnames(variance) <-
    do.call("paste", c(lapply(1:ncol(valuesout), function(i) {
      paste(names(valuesout)[i], "=", round(valuesout[[i]], 2))
    }), sep = "; "))

  valuesout$est <- est
  valuesout$se <- sqrt(diag(variance))


  structure(list(estimates = valuesout, covariance = variance),
            class = "std")

}
