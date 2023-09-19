#' @title Get regression standardized estimates from a glm
#' @param formula
#' The formula which is used to fit the glm model for the outcome.
#' @param data The data.
#' @param family
#' The family argument which is used to fit the glm model for the outcome.
#' @param values
#' A named list or data.frame specifying the variables and values
#' at which marginal means of the outcome will be estimated.
#' @param clusterid
#' An optional string containing the name of a cluster identification variable
#' when data are clustered.
#' @param case_control
#' Whether the data comes from a case-control study.
#' @param p_population
#' Specifies the incidence in the population when \code{case_control=TRUE}.
#' @param matched_density_cases
#' A function of the matching variable.
#' The probability (or density) of the matched variable among the cases.
#' @param matched_density_controls
#' A function of the matching variable.
#' The probability (or density) of the matched variable among the controls.
#' @param matching_variable
#' The matching variable extracted from the data set.
#' @param ci_type A string, indicating the type of confidence intervals.
#' Either "plain", which gives untransformed intervals, or "log", which gives
#' log-transformed intervals.
#' @param ci_level Coverage probability of confidence intervals.
#' @param transforms A vector of transforms in the following format:
#' If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized
#' mean \eqn{\theta(x)} is transformed into \eqn{\psi(x)=log\{\theta(x)\}},
#' \eqn{\psi(x)=log[\theta(x)/\{1-\theta(x)\}]}, or
#' \eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively.
#' If the vector is \code{NULL}, then \eqn{\psi(x)=\theta(x)}.
#' @param contrasts A vector of contrasts in the following format:
#' If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(x)-\psi(x_0)}
#' or \eqn{\psi(x) / \psi(x_0)} are constructed, where \eqn{x_0} is a reference
#' level specified by the \code{reference} argument. Has to be be \code{NULL}
if no references are specified.
#' @param references A vector of references in the following format:
#' If \code{contrasts} is not \code{NULL}, the desired reference level(s).
#' @returns
#' An object of class \code{std_glm}.
#' This is basically a list with components estimates and covariance estimates.
#' The output contains estimates for contrasts and confidence intervals
#' for all combinations of \code{transforms}, \code{references}
# 'and \code{transforms}.
#' @details \code{standardize_glm} performs regression standardization
#' in generalized linear models,
#' at specified values of the exposure, over the sample covariate distribution.
#' Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome, the exposure, and a
#' vector of covariates, respectively.
#' \code{standardize_glm} uses a fitted generalized linear
#' model to estimate the standardized
#' mean \eqn{\theta(x)=E\{E(Y|X=x,Z)\}},
#' where \eqn{x} is a specific value of \eqn{X},
#' and the outer expectation is over the marginal distribution of \eqn{Z}.
#' @references Rothman K.J., Greenland S., Lash T.L. (2008).
#' \emph{Modern Epidemiology}, 3rd edition.
#' Lippincott, Williams \& Wilkins.
#' @references Sjolander A. (2016).
#' Regression standardization with the R-package stdReg.
#' \emph{European Journal of Epidemiology} \bold{31}(6), 563-574.
#' @references Sjolander A. (2016).
#' Estimation of causal effect measures with the R-package stdReg.
#' \emph{European Journal of Epidemiology} \bold{33}(9), 847-858.
#' @examples
#'
#' # basic example
#' # needs to correctly specify the outcome model and no unmeasered confounders
#' # (+ standard causal assunmptions)
#' set.seed(6)
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
#' dd <- data.frame(Z, X, Y)
#' x <- standardize_glm(
#'   formula = Y ~ X * Z,
#'   family = "binomial",
#'   data = dd,
#'   values = list(X = 0:1),
#'   contrasts = c("difference", "ratio"),
#'   reference = 0
#' )
#' x
#' # different transformations of causal effects
#'
#' # example from Sjolander (2016) with case-control data
#' # here the matching variable needs to be passed as an argument
#' library(AF)
#' data("singapore")
#' Mi <- singapore$Age
#' m <- mean(Mi)
#' s <- sd(Mi)
#' d <- 5
#' standardize_glm(
#'   formula = Oesophagealcancer ~ (Everhotbev + Age + Dial + Samsu + Cigs)^2,
#'   family = binomial, data = singapore,
#'   values = list(Everhotbev = 0:1), clusterid = "Set",
#'   case_control = TRUE,
#'   matched_density_cases = function(x) dnorm(x, m, s),
#'   matched_density_controls = function(x) dnorm(x, m - d, s),
#'   matching_variable = Mi,
#'   p_population = 19.3 / 100000
#' )
#'
#' # multiple exposures
#' set.seed(7)
#' n <- 100
#' Z <- rnorm(n)
#' X1 <- rnorm(n, mean = Z)
#' X2 <- rnorm(n)
#' Y <- rbinom(n, 1, prob = (1 + exp(X1 + X2 + Z))^(-1))
#' dd <- data.frame(Z, X1, X2, Y)
#' x <- standardize_glm(
#'   formula = Y ~ X1 + X2 + Z,
#'   family = "binomial",
#'   data = dd, values = list(X1 = 0:1, X2 = 0:1),
#'   contrasts = c("difference", "ratio"),
#'   references = "0, 0"
#' )
#' x
#'
#' # continuous exposure
#' set.seed(2)
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' Y <- rnorm(n, mean = X + Z + 0.1 * X^2)
#' dd <- data.frame(Z, X, Y)
#' x <- standardize_glm(
#'   formula = Y ~ X * Z,
#'   family = "gaussian",
#'   data = dd,
#'   values = list(X = seq(-1, 1, 0.1))
#' )
#'
#' # plot standardized mean as a function of x
#' plot(x)
#' # plot standardized mean - standardized mean at x = 0 as a function of x
#' plot(x, contrast = "difference", reference = 0)
#'
#' @export standardize_glm
standardize_glm <- function(formula,
                            data,
                            values,
                            clusterid,
                            matched_density_cases,
                            matched_density_controls,
                            matching_variable,
                            p_population,
                            case_control = FALSE,
                            ci_level = 0.95,
                            ci_type = "plain",
                            contrasts = NULL,
                            family = "gaussian",
                            references = NULL,
                            transforms = NULL) {
  n <- nrow(data)

  if (!inherits(values, c("data.frame", "list"))) {
    stop("values is not an object of class list or data.frame")
  }

  ## Check that the names of values appear in the data
  check_values_data(values, data)

  ## Set various relevant variables
  valuesout <- outcome <- exposure_names <- NULL
  list2env(get_outcome_exposure(formula, data, values), envir = environment())

  if (case_control) {
    if (missing(p_population)) {
      stop("you have to specify population prevalence when case_control = TRUE")
    }
    if (!(identical(family, binomial)) || !is.binary(outcome)) {
      stop("the option 'case_control = TRUE' only works with a binary outcome")
    }
    cases <- which(outcome == 1L)
    controls <- which(outcome == 0L)
    n1 <- length(cases)
    p_star <- n1 / n
    n0 <- n - n1
    weights <- outcome * p_population / p_star +
      (1.0 - outcome) * (1.0 - p_population) / (1.0 - p_star) *
        matched_density_controls(matching_variable) /
        matched_density_cases(matching_variable)
  } else {
    weights <- rep(1.0, n)
  }
  data[["weights"]] <- weights
  fit_outcome <- fit_glm(formula, family, data, "outcome")

  ## Estimation and variance estimation
  ## In the implementation for the Score equations,
  ## the Score equations corresponding to the standardized mean come first,
  ## with the the Score equations corresponding
  ## for the parameters of the glm coming afterwards

  ## Get the derivative of the inverse link function
  deriv_inv_link <- fit_outcome[["family"]][["mu.eta"]]
  ## Contains predictions for data where the values of the data has been replaced by the i'th row of valuesout
  predmat <- matrix(nrow = n, ncol = nrow(valuesout))
  ## Is used for the variance calculation where it corresponds to the derivative of the EE of theta(x) wrt. beta
  dmu_dbeta <- matrix(nrow = nrow(valuesout), ncol = length(fit_outcome[["coefficients"]]))
  for (i in seq_len(nrow(valuesout))) {
    ## Get the data for use with the predictions,
    ## i.e., the data from the fit but the with the variables in values replaced
    ## with the i'th row of valuesout
    data_x <- do.call("transform", c(
      list(data),
      valuesout[i, , drop = FALSE]
    ))

    ## Save the predictions for data_x
    pred_x <- predict(object = fit_outcome, newdata = data_x)
    predmat[, i] <- fit_outcome[["family"]][["linkinv"]](pred_x)

    ## Corresponds to the derivative of eta^{-1}(h(X,Z; beta)) wrt. beta, i.e., use chain rule
    dmu_dbeta[i, ] <- colMeans(weights * deriv_inv_link(pred_x) * model.matrix(object = terms(fit_outcome), data = data_x))
  }
  ## Estimates of standardized means
  estimates <- colSums(weights * predmat, na.rm = TRUE) / sum(weights)

  ## Implement variance estimation according to Appendix 1 of Sjolander, A. (2016)

  ## Get Score (U) and the Fisher information matrix (I) from the glm object
  ## NOTE: The dispersion parameter is set to 1 for the variance calculations
  sandwich_fit <- sandwich(fit = fit_outcome, data = data, weights = weights)

  ## Estimate the term that "corresponds" to var{U_{v,i}(\nu)} of Equation (5)
  ee_beta <- sandwich_fit[["U"]]
  ee_means <- weights * (predmat - matrix(rep(estimates, each = n),
    nrow = n,
    ncol = nrow(valuesout)
  )) ## EE corresponding to the standardized means (or rather each term in estimating equation)
  ee <- cbind(ee_means, ee_beta)
  if (missing(clusterid)) {
    if (case_control) {
      j_mat <- n0 / n * var(ee[controls, ], na.rm = TRUE) + n1 / n * var(ee[cases, ], na.rm = TRUE)
    } else {
      j_mat <- var(ee, na.rm = TRUE)
    }
  } else {
    clusters <- data[, clusterid]
    ee <- aggr(x = ee, clusters = clusters)
    if (case_control) {
      warning("case_control = TRUE may not give reasonable results for the variance with clustering")
      ## Maybe we need adjust the variance as we do above with no clustering?
    }
    j_mat <- var(ee, na.rm = TRUE)
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

  upper_i_mat <- cbind(-diag(nrow(valuesout)) * mean(weights), dmu_dbeta)
  lower_i_mat <- cbind(matrix(0.0,
    nrow = length(fit_outcome[["coefficients"]]),
    ncol = nrow(valuesout)
  ), sandwich_fit[["I"]])
  i_mat <- rbind(upper_i_mat, lower_i_mat)

  ## Apply Equation (5) of Sjolander, A. (2016)
  if (missing(clusterid)) {
    variance <- (solve(i_mat) %*% j_mat %*% t(solve(i_mat)) / n)[seq_len(nrow(valuesout)), seq_len(nrow(valuesout))]
  } else {
    variance <- (solve(i_mat) %*% j_mat %*% t(solve(i_mat)) * length(unique(clusters)) / n^2L)[seq_len(nrow(valuesout)), seq_len(nrow(valuesout))]
  }
  fit_exposure <- NULL

  format_result_standardize_glm(
    contrasts,
    references,
    transforms,
    valuesout,
    variance,
    fit_outcome,
    fit_exposure,
    exposure_names,
    estimates,
    ci_type,
    ci_level
  )
}

#' @title Get regression standardized doubly-robust estimates from a glm
#' @param formula_outcome
#' The formula which is used to fit the glm model for the outcome.
#' @param formula_exposure
#' The formula which is used to fit the glm model for the exposure.
#' If not \code{NULL},
#' a doubly robust estimator of the standardized estimator is used.
#' @param family_outcome
#' The family argument which is used to fit the glm model for the outcome.
#' @param family_exposure
#' The family argument which is used to fit the glm model for the exposure.
#' @inherit standardize_glm
#' @details \code{standardize_glm_dr} performs regression standardization
#' in generalized linear models, see e.g., documentation for \code{standardize_glm_dr}. Specifically,
#' this version uses a doubly robust estimator for standardization, meaning inference is valid
#' when either the outcome regression or the exposure model is correctly specified
#' and there is no unmeasered confounding.
#' @references Gabriel E.E., Sachs, M.C., Martinussen T., Waernbaum I.,
#' Goetghebeur E., Vansteelandt S., Sjolander A. (????),
#' Inverse probability of treatment weighting with
#' generalized linear outcome models for doubly robust estimation.
#' ????
#' @examples
#'
#' # doubly robust estimator
#' # needs to correctly specify either the outcome model or the exposure model
#' # for confounding
#' # NOTE: only works with binary exposures
#' library(AF)
#' data <- AF::clslowbwt
#' x <- standardize_glm_dr(
#'   formula_outcome = bwt ~ smoker * (race + age + lwt) + I(age^2) + I(lwt^2),
#'   formula_exposure = smoker ~ race * age * lwt + I(age^2) + I(lwt^2),
#'   family_outcome = "gaussian",
#'   family_exposure = "binomial",
#'   data = data,
#'   values = list(smoker = c(0, 1)), contrasts = "difference", reference = 0
#' )
#'
#' set.seed(6)
#' n <- 100
#' Z <- rnorm(n)
#' X <- rbinom(n, 1, prob = (1 + exp(Z))^(-1))
#' Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
#' dd <- data.frame(Z, X, Y)
#' x <- standardize_glm_dr(
#'   formula_outcome = Y ~ X * Z, formula_exposure = X ~ Z,
#'   family_outcome = "binomial",
#'   data = dd,
#'   values = list(X = 0:1), references = c(0, 1), contrasts = c("difference"), transforms = c("odds")
#' )
#'
#' @export standardize_glm_dr
standardize_glm_dr <- function(formula_exposure,
                               formula_outcome,
                               data,
                               values,
                               ci_level = 0.95,
                               ci_type = "plain",
                               contrasts = NULL,
                               family_outcome = "gaussian",
                               family_exposure = "binomial",
                               references = NULL,
                               transforms = NULL) {
  ## Preparation and various checks
  n <- nrow(data)

  if (!inherits(values, c("data.frame", "list"))) {
    stop("values is not an object of class list or data.frame")
  }

  ## Check that the names of values appear in the data
  check_values_data(values, data)

  ## Set various relevant variables
  valuesout <- outcome <- exposure_names <- exposure <- NULL
  list2env(get_outcome_exposure(formula_outcome, data, values), envir = environment())

  if (length(exposure_names) > 1L) {
    stop("there has to be only one exposure
              with the doubly robust estimator")
  }

  ## Check that exposure is binary
  if (inherits(family_exposure, "function") &&
    !(identical(family_exposure, binomial)) ||
    (inherits(family_exposure, "character") && family_exposure != "binomial") ||
    !is.binary(exposure) || nrow(valuesout) != 2L) {
    stop("the exposure has to be binary (0 or 1)")
  }
  data[["weights"]] <- rep(1, n)
  fit_exposure <- fit_glm(formula_exposure, family_exposure, data, "exposure")

  g_weights <- predict(fit_exposure, type = "response")
  data[["weights"]] <- exposure / g_weights + (1.0 - exposure) / (1.0 - g_weights)
  fit_outcome <- fit_glm(formula_outcome, family_outcome, data, "outcome")

  ## Create two data sets with every observation treated or untreated
  data_exposure_1 <- data_exposure_0 <- data
  data_exposure_1[[exposure_names]] <- 1L
  data_exposure_0[[exposure_names]] <- 0L

  ## Calculate estimates
  standardized_estimate_1 <- mean(predict(fit_outcome, newdata = data_exposure_1, type = "response"))
  standardized_estimate_0 <- mean(predict(fit_outcome, newdata = data_exposure_0, type = "response"))

  ## Variance estimation based on Appendix of ????

  ## refit using unweighted estimating equations for variance estimation
  fit_outcome_unweighted <- glm(formula_outcome, family = family_outcome, data = data)
  est1 <- predict(fit_outcome_unweighted, newdata = data_exposure_1, type = "response")
  est0 <- predict(fit_outcome_unweighted, newdata = data_exposure_0, type = "response")
  mme <- model.matrix(fit_exposure)
  mmoe <- mmou <- model.matrix(fit_outcome)
  mmoe[, exposure_names] <- 1L
  mmou[, exposure_names] <- 0L

  ## IF for the parametric models
  if_exposure <- t(vcov(fit_exposure) %*% t(sandwich::estfun(fit_exposure)))
  if_outcome <- t(vcov(fit_outcome) %*% t(sandwich::estfun(fit_outcome_unweighted)))

  eif_terms_1 <- (exposure / g_weights * (outcome - est1) +
    (est1 - standardized_estimate_1)) / n
  eif_terms_0 <- ((1.0 - exposure) / (1.0 - g_weights) * (outcome - est0) +
    (est0 - standardized_estimate_0)) / n
  g_dot <- family(fit_exposure)[["mu.eta"]](predict(fit_exposure, type = "link"))
  r_dot_0 <- family(fit_outcome)[["mu.eta"]](predict(fit_outcome_unweighted, newdata = data_exposure_0, type = "link"))
  r_dot_1 <- family(fit_outcome)[["mu.eta"]](predict(fit_outcome_unweighted, newdata = data_exposure_1, type = "link"))

  ## NOTE: chain rule is used for gi.dot and ri.dot and signs for terms corresponding to unexposed are changed
  k_term_1 <- (-1.0 / n) * matrix(((exposure * g_dot) / g_weights^2L) *
    (outcome - est1), nrow = 1L, ncol = n) %*% mme
  k_term_0 <- (1.0 / n) * matrix((((1.0 - exposure) * g_dot) / (1.0 - g_weights)^2L) *
    (outcome - est0), nrow = 1L, ncol = n) %*% mme

  ## NOTE: we clearly have to use the unweighted outcome fits, so that the Lterms correspond with the one in the article
  l_term_1 <- (1.0 / n) * matrix(r_dot_1 * (1.0 - exposure / g_weights),
    nrow = 1L, ncol = n
  ) %*% mmoe
  l_term_0 <- (1.0 / n) * matrix(r_dot_0 * ((1.0 - exposure) / (1.0 - g_weights) - 1.0),
    nrow = 1L, ncol = n
  ) %*% mmou

  if_1 <- rowSums(cbind(eif_terms_1, (if_exposure %*% t(k_term_1)), (if_outcome %*% t(l_term_1))))
  if_0 <- rowSums(cbind(eif_terms_0, (if_exposure %*% t(k_term_0)), (if_outcome %*% t(l_term_0))))
  covar <- sum(if_0 * if_1)
  variance <- matrix(c(sum(if_0^2L), covar, covar, sum(if_1^2L)), nrow = 2L)
  estimates <- c(standardized_estimate_0, standardized_estimate_1)

  format_result_standardize_glm(
    contrasts,
    references,
    transforms,
    valuesout,
    variance,
    fit_outcome,
    fit_exposure,
    exposure_names,
    estimates,
    ci_type,
    ci_level
  )
}

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

fit_glm <- function(formula, family, data, response) {
  weights <- NULL
  ## fit with quasipoisson/quasibinomial to suppress warnings
  if (response == "outcome") {
    if ((inherits(family, "function") && identical(family, stats::binomial)) ||
      (inherits(family, "character") && family == "binomial")) {
      family <- stats::quasibinomial
    } else if ((inherits(family, "function") && identical(family, stats::poisson)) ||
      (inherits(family, "character") && family == "poisson")) {
      family <- stats::quasipoisson
    }
  }

  ## try fitting a glm model
  fit <- tryCatch(
    {
      glm(formula = formula, data = data, family = family, weights = weights)
    },
    error = function(cond) {
      return(cond)
    }
  )
  if (inherits(fit, "simpleError")) {
    stop("glm for ", response, " failed with error: ", fit[["message"]])
  } else {
    fit
  }
}

format_result_standardize_glm <- function(contrasts,
                                          references,
                                          transforms,
                                          valuesout,
                                          variance,
                                          fit_outcome,
                                          fit_exposure,
                                          exposure_names,
                                          estimates,
                                          ci_type,
                                          ci_level) {
  contrast <- reference <- NULL
  ## Add names to asymptotic covariance matrix
  rownames(variance) <- colnames(variance) <-
    do.call("paste", c(lapply(seq_len(ncol(valuesout)), function(i) {
      paste(names(valuesout)[i], "=", round(valuesout[[i]], 2L))
    }), sep = "; "))

  valuesout[["se"]] <- sqrt(diag(variance))
  valuesout[["estimates"]] <- estimates
  res <- structure(
    list(
      estimates = valuesout,
      covariance = variance,
      fit_outcome = fit_outcome,
      fit_exposure = fit_exposure,
      exposure_names = exposure_names
    ),
    class = "std_glm_helper"
  )
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
    summary.std_glm_helper(res,
      ci_type = ci_type,
      ci_level = ci_level,
      transform = transform,
      contrast = contrast,
      reference = reference
    )
  }
  res_contrast <- as.list(as.data.frame(do.call(mapply, c("summary_fun", unname(as.list(grid))))))
  res_fin <- list(res_contrast = res_contrast, res = res)
  class(res_fin) <- "std_glm"
  res_fin
}

summary.std_glm_helper <- function(object, ci_type = "plain", ci_level = 0.95,
                                   transform = NULL, contrast = NULL, reference = NULL, ...) {
  null_helper <- function(x) {
    if (is.null(x) || x == "NULL") {
      NULL
    } else {
      x
    }
  }
  transform <- null_helper(transform)
  contrast <- null_helper(contrast)
  reference <- null_helper(reference)

  est_old_table <- object[["estimates"]]
  est <- est_old_table[["estimates"]]
  v_mat <- as.matrix(object[["covariance"]])
  n_x_levs <- nrow(est_old_table)
  if (!is.null(transform)) {
    if (transform == "log") {
      if (any(est <= 0)) {
        stop("transform='log' requires that the (standardized) estiamtes are positive.")
      }
      dtransform_dm <- diag(1.0 / est, nrow = n_x_levs, ncol = n_x_levs)
      est <- log(est)
    }
    if (transform == "logit") {
      if (any(est <= 0 | est >= 1)) {
        stop("transform='logit' requires that the (standardized) estimates take values in (0, 1).")
      }
      dtransform_dm <- diag(1.0 / (est * (1.0 - est)), nrow = n_x_levs, ncol = n_x_levs)
      est <- logit(est)
    }
    if (transform == "odds") {
      if (any(est == 1)) {
        stop("transform='odds' requires that the (standardized) estimates are not equal to 1. ")
      }
      dtransform_dm <- diag(1.0 / (1.0 - est)^2L, nrow = n_x_levs, ncol = n_x_levs)
      est <- odds(est)
    }
    v_mat <- t(dtransform_dm) %*% v_mat %*% dtransform_dm
  }
  if (!is.null(contrast)) {
    if (is.null(reference)) {
      stop("When specifying contrast, reference must be specified as well")
    }
    reference <- gsub(" ", "", reference, fixed = TRUE)
    if (length(object[["exposure_names"]]) > 1L) {
      est_old_table[["exposure"]] <- do.call(paste, c(est_old_table[, object[["exposure_names"]]], sep = ","))
    } else {
      est_old_table[["exposure"]] <- est_old_table[[object[["exposure_names"]]]]
    }

    referencepos <- match(reference, est_old_table[["exposure"]])
    if (is.na(referencepos)) {
      stop("reference must be a value in x")
    }
    if (contrast == "difference") {
      dcontrast_dtransform <- diag(n_x_levs)
      dcontrast_dtransform[referencepos, ] <- -1.0
      dcontrast_dtransform[referencepos, referencepos] <- 0.0
      est <- est - est[referencepos]
    } else if (contrast == "ratio") {
      dcontrast_dtransform <- diag(1.0 / est[referencepos], nrow = n_x_levs, ncol = n_x_levs)
      dcontrast_dtransform[referencepos, ] <- -est / est[referencepos]^2L
      dcontrast_dtransform[referencepos, referencepos] <- 1.0
      est <- est / est[referencepos]
    } else {
      stop("contrast not supported.")
    }
    v_mat <- t(dcontrast_dtransform) %*% v_mat %*% dcontrast_dtransform
    v_mat[referencepos, ] <- 0.0
    v_mat[, referencepos] <- 0.0
  }

  var <- diag(v_mat)
  se <- sqrt(var)
  conf_int <- CI(est = est, var = var, CI.type = ci_type, CI.level = ci_level)

  if (is.factor(reference)) {
    reference <- as.character(reference)
  }
  if (!is.null(contrast)) {
    est_table <- data.frame(est_old_table[["exposure"]], as.matrix(cbind(est, se, conf_int), nrow = length(est), ncol = 4L))
    colnames(est_table) <- c("Exposure", "Estimate", "Std. Error", paste("lower", ci_level), paste("upper", ci_level))
  } else {
    est_table <- data.frame(est_old_table[, object[["exposure_names"]]], as.matrix(cbind(est, se, conf_int), nrow = length(est), ncol = 4L))
    colnames(est_table) <- c(object[["exposure_names"]], "Estimate", "Std. Error", paste("lower", ci_level), paste("upper", ci_level))
  }
  rownames(est_table) <- NULL
  out <- c(object, list(
    est_table = est_table, transform = transform,
    contrast = contrast, reference = reference,
    ci_type = ci_type, ci_level = ci_level
  ))

  class(out) <- "summary.std_glm_helper"
  return(out)
}

#' @title Prints summary of GLM regression standardization fit
#' @param x an object of class \code{"std_glm"}.
#' @param \dots unused
#' @rdname print
#' @export print.std_glm
#' @export
print.std_glm <- function(x, ...) {
  if (!is.null(x[["res"]][["fit_exposure"]])) {
    cat("Doubly robust estimator with: \n")
    cat("\nExposure formula: ")
    print(x[["res"]][["fit_exposure"]][["formula"]])
    cat("Exposure link function:", x[["res"]][["fit_exposure"]][["family"]][["link"]], "\n")
  }
  cat("Outcome formula: ")
  print(x[["res"]][["fit_outcome"]][["formula"]])
  cat("Outcome family:", x[["res"]][["fit_outcome"]][["family"]][["family"]], "\n")
  cat("Outcome link function:", x[["res"]][["fit_outcome"]][["family"]][["link"]], "\n")
  cat("Exposure: ", toString(x[["res"]][["exposure_names"]]), "\n")
  cat("\n")
  cat("Tables: \n")
  for (l in seq_len(length(x[["res_contrast"]]))) {
    temp <- x[["res_contrast"]][[paste0("V", l)]]
    if (!is.null(temp[["transform"]])) {
      cat("Transform: ", levels(temp[["transform"]])[[temp[["transform"]]]], "\n")
    }
    if (!is.null(temp[["contrast"]])) {
      cat("Reference level: ", temp[["input"]][["X"]], "=", temp[["reference"]], "\n")
      cat("Contrast: ", levels(temp[["contrast"]])[[temp[["contrast"]]]], "\n")
    }
    print(temp[["est_table"]], digits = 3L)
    cat("\n")
  }
}

#' @title Plots GLM regression standardization fit
#' @description This is a \code{plot} method for class \code{"std_glm"}.
#' @param x An object of class \code{"std_glm"}
#' @param ci_type A string, indicating the type of confidence intervals. Either "plain", which
#' gives untransformed intervals, or "log", which gives log-transformed intervals.
#' @param ci_level Coverage probability of confidence intervals.
#' @param transform  If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized
#' mean \eqn{\theta(x)} is transformed into \eqn{\psi(x)=log\{\theta(x)\}},
#' \eqn{\psi(x)=log[\theta(x)/\{1-\theta(x)\}]}, or
#' \eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively. If left unspecified,
#' \eqn{\psi(x)=\theta(x)}.
#' @param contrast If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(x)-\psi(x_0)}
#' or \eqn{\psi(x) / \psi(x_0)} are constructed, where \eqn{x_0} is a reference
#' level specified by the \code{reference} argument.
#' If not \code{NULL}, a doubly robust estimator of the standardized estimator is used.
#' @param reference If \code{contrast} is specified, the desired reference level.
#' @param \dots Unused.
#' @examples
#' # see standardize_glm
#'
#' @rdname plot
#' @export plot.std_glm
#' @export
plot.std_glm <- function(x, ci_type = "plain", ci_level = 0.95,
                         transform = NULL, contrast = NULL, reference = NULL, ...) {
  object <- x[["res"]]
  x <- object[["estimates"]][, object[["exposure_names"]]]

  dots <- list(...)

  xlab <- object[["exposure_names"]]

  if (length(xlab) > 1L) {
    stop("cannot do plot with multiple exposures")
  }

  if (is.factor(reference)) {
    reference <- as.character(reference)
  }

  if (is.null(contrast)) {
    if (is.null(transform)) {
      ylab <- expression(mu)
    } else {
      if (transform == "log") {
        ylab <- expression(paste("log(", mu, ")"))
      }
      if (transform == "logit") {
        ylab <- expression(paste("logit(", mu, ")"))
      }
      if (transform == "odds") {
        ylab <- expression(paste(mu, "/(1-", mu, ")"))
      }
    }
  } else {
    if (contrast == "difference") {
      if (is.null(transform)) {
        ylab <- c(bquote(paste(mu, "-", mu[.(reference)])), expression())
      } else {
        if (transform == "log") {
          ylab <- c(bquote(paste(
            log, "(", mu, ")-", log, "(",
            mu[.(reference)], ")"
          )), expression())
        }
        if (transform == "logit") {
          ylab <- c(bquote(paste(
            logit, "(", mu, ")-", logit,
            "(", mu[.(reference)], ")"
          )), expression())
        }
        if (transform == "odds") {
          ylab <- c(
            bquote(paste(
              mu, "/(1-", mu, ")-",
              mu[.(reference)], "/(1-", mu[.(reference)], ")"
            )),
            expression()
          )
        }
      }
    }
    if (contrast == "ratio") {
      if (is.null(transform)) {
        ylab <- c(bquote(paste(mu, "/", mu[.(reference)])), expression())
      } else {
        if (transform == "log") {
          ylab <- c(bquote(paste(
            log, "(", mu, ")/", log, "(",
            mu[.(reference)], ")"
          )), expression())
        }
        if (transform == "logit") {
          ylab <- c(bquote(paste(
            logit, "(", mu, ")/", logit,
            "(", mu[.(reference)], ")"
          )), expression())
        }
        if (transform == "odds") {
          ylab <- c(bquote(paste(
            mu, "/(1-", mu, ")/",
            mu[.(reference)], "/(1-", mu[.(reference)],
            ")"
          )), expression())
        }
      }
    }
  }

  sum.obj <- summary(
    object = object, ci_type = ci_type, ci_level = ci_level,
    transform = transform, contrast = contrast, reference = reference
  )
  est <- sum.obj[["est_table"]][, 2L]
  lower <- sum.obj[["est_table"]][, 4L]
  upper <- sum.obj[["est_table"]][, 5L]
  ylim <- c(min(c(lower, upper)), max(c(lower, upper)))
  if (is.numeric(x) && length(x) > 1L) {
    args <- list(x = x, y = x, xlab = xlab, ylab = ylab, ylim = ylim, type = "n")
    args[names(dots)] <- dots
    do.call("plot", args = args)
    lines(x, est)
    lines(x, upper, lty = 3L)
    lines(x, lower, lty = 3L)
  }
  if (is.factor(x) || is.binary(x) || (is.numeric(x) && length(x) == 1L)) {
    x_seq <- seq_len(length(x))
    args <- list(
      x = x_seq, y = seq_len(length(x)), xlab = xlab, ylab = ylab,
      xlim = c(0.0, length(x) + 1.0), ylim = ylim, type = "n", xaxt = "n"
    )
    args[names(dots)] <- dots
    do.call("plot", args = args)
    points(x_seq, est)
    points(x_seq, upper, pch = 0L)
    points(x_seq, lower, pch = 0L)
    for (i in x_seq) {
      lines(x = c(i, i), y = c(lower[[i]], upper[[i]]), lty = "dashed")
    }
    mtext(text = x, side = 1L, at = x_seq)
  }
}
