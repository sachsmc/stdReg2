#' @title Get standardized estimates using the g-formula with a custom model
#' @inherit standardize_glm
#' @param fitter The function to call to fit the data.
#' @param arguments
#' The arguments to be used in the fitter function as a \code{list}.
#' @param predict_fun The function used to predict the means/probabilities
#' for a new data set on the response level. For survival data,
#' this should be a matrix where each column is the time, and each
#' row the data.
#' @param times For use with survival data. Set to \code{NULL} otherwise.
#' @param B Number of nonparametric bootstrap resamples. Default is \code{NULL} (no bootstrap).
#' @param seed The seed to use with the nonparametric bootstrap.
#' @param progressbar Logical, if TRUE will print bootstrapping progress to the console
#' @returns
#' An object of class \code{std_custom}.
#' This is a list with components estimates and fit for the outcome model.
#' @details
#' Let \eqn{Y}, \eqn{X}, and \eqn{Z} be the outcome, the exposure, and a
#' vector of covariates, respectively.
#' \code{standardize} uses a
#' model to estimate the standardized
#' mean \eqn{\theta(x)=E\{E(Y|X=x,Z)\}},
#' where \eqn{x} is a specific value of \eqn{X},
#' and the outer expectation is over the marginal distribution of \eqn{Z}.
#' With survival data, \eqn{Y=I(T > t)},
#' and a vector of different time points \code{times} (\eqn{t}) can be given,
#' where \eqn{T} is the uncensored survival time.
#' @examples
#'
#' set.seed(6)
#' n <- 100
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' Y <- rbinom(n, 1, prob = (1 + exp(X + Z))^(-1))
#' dd <- data.frame(Z, X, Y)
#' prob_predict.glm <- function(...) predict.glm(..., type = "response")
#'
#' x <- standardize(
#'   fitter = "glm",
#'   arguments = list(
#'     formula = Y ~ X * Z,
#'     family = "binomial"
#'   ),
#'   predict_fun = prob_predict.glm,
#'   data = dd,
#'   values = list(X = seq(-1, 1, 0.1)),
#'   B = 100,
#'   reference = 0,
#'   contrasts = "difference"
#' )
#' x
#'
#' require(survival)
#' prob_predict.coxph <- function(object, newdata, times) {
#'   fit.detail <- suppressWarnings(basehaz(object))
#'   cum.haz <- fit.detail$hazard[sapply(times, function(x) max(which(fit.detail$time <= x)))]
#'   predX <- predict(object = object, newdata = newdata, type = "risk")
#'   res <- matrix(NA, ncol = length(times), nrow = length(predX))
#'   for (ti in seq_len(length(times))) {
#'     res[, ti] <- exp(-predX * cum.haz[ti])
#'   }
#'   res
#' }
#' set.seed(68)
#' n <- 500
#' Z <- rnorm(n)
#' X <- rnorm(n, mean = Z)
#' T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
#' C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
#' U <- pmin(T, C) # time at risk
#' D <- as.numeric(T < C) # event indicator
#' dd <- data.frame(Z, X, U, D)
#' x <- standardize(
#' fitter = "coxph",
#'   arguments = list(
#'     formula = Surv(U, D) ~ X + Z + X * Z,
#'     method = "breslow",
#'     x = TRUE,
#'     y = TRUE
#'   ),
#'   predict_fun = prob_predict.coxph,
#'   data = dd,
#'   times = 1:5,
#'   values = list(X = c(-1, 0, 1)),
#'   B = 100,
#'   reference = 0,
#'   contrasts = "difference"
#' )
#' x
#' @export standardize
standardize <- function(fitter,
                        arguments,
                        predict_fun,
                        data,
                        values,
                        B = NULL,
                        ci_level = 0.95,
                        contrasts = NULL,
                        reference = NULL,
                        seed = NULL,
                        times = NULL,
                        transforms = NULL,
                        progressbar = TRUE) {
  ## Preparation and various checks
  n <- nrow(data)

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
  exposure_names <- colnames(valuesout)

  fit_outcome <- fit_helper(arguments, fitter, data)
  estimate_fun <- function(valuesout, times, data, fit_outcome) {
    if (is.null(times)) {
      estimates <- rep(NA, nrow(valuesout))
    } else {
      estimates <- matrix(NA, ncol = length(times), nrow = nrow(valuesout))
      colnames(estimates) <- paste("t =", times)
    }

    for (i in seq_len(nrow(valuesout))) {
      data_x <- do.call("transform", c(
        list(data),
        valuesout[i, , drop = FALSE]
      ))

      ## Save the predictions for data_x
      if (is.null(times)) {
        estimates[i] <- mean(predict_fun(object = fit_outcome, newdata = data_x))
      } else {
        estimates[i, ] <- colMeans(predict_fun(object = fit_outcome, newdata = data_x, times = times))
      }
    }
    estimates
  }
  estimates <- estimate_fun(valuesout, times, data, fit_outcome)
  estimates_boot <- list()
  if (!is.null(B)) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    pb <- utils::txtProgressBar(
      min = 1,
      max = B,
      style = 3,
      width = 50
    )

    #cat("Bootstrapping... This may take some time... \n")
    for (b in seq_len(B)) {
     if(progressbar) utils::setTxtProgressBar(pb, b)
      data_boot <- data[sample(seq_len(n), replace = TRUE), ]
      fit_outcome_boot <- fit_helper(arguments, fitter, data_boot)
      estimates_boot[[b]] <- estimate_fun(valuesout, times, data_boot, fit_outcome_boot)
    }
  }
  valuesout <- cbind(valuesout, estimates)
  res <- structure(
    list(
      B = B,
      estimates = valuesout,
      fit_outcome = fit_outcome,
      estimates_boot = estimates_boot,
      exposure_names = exposure_names,
      times = times
    )
  )
  format_result_standardize(
    res,
    contrasts,
    reference,
    transforms,
    "plain",
    ci_level,
    "std_custom",
    "summary_standardize"
  )
}

#' @title Get standardized estimates using the g-formula with and separate models for each exposure level in the data
#' @inherit standardize
#' @param fitter_list The function to call to fit the data (as a list).
#' @param predict_fun_list The function used to predict the means/probabilities
#' for a new data set on the response level. For survival data,
#' this should be a matrix where each column is the time, and each
#' row the data (as a list).
#' @details
#' See \code{standardize}. The difference is here that different models
#' can be fitted for each value of \code{x} in \code{values}.
#' @examples
#'
#' require(survival)
#' prob_predict.coxph <- function(object, newdata, times) {
#'   fit.detail <- suppressWarnings(basehaz(object))
#'   cum.haz <- fit.detail$hazard[sapply(times, function(x) max(which(fit.detail$time <= x)))]
#'   predX <- predict(object = object, newdata = newdata, type = "risk")
#'   res <- matrix(NA, ncol = length(times), nrow = length(predX))
#'   for (ti in seq_len(length(times))) {
#'     res[, ti] <- exp(-predX * cum.haz[ti])
#'   }
#'   res
#' }
#'
#' set.seed(68)
#' n <- 500
#' Z <- rnorm(n)
#' X <- rbinom(n, 1, prob = 0.5)
#' T <- rexp(n, rate = exp(X + Z + X * Z)) # survival time
#' C <- rexp(n, rate = exp(X + Z + X * Z)) # censoring time
#' U <- pmin(T, C) # time at risk
#' D <- as.numeric(T < C) # event indicator
#' dd <- data.frame(Z, X, U, D)
#' x <- standardize_level(
#'   fitter_list = list("coxph", "coxph"),
#'   arguments = list(
#'     list(
#'       formula = Surv(U, D) ~ X + Z + X * Z,
#'       method = "breslow",
#'       x = TRUE,
#'       y = TRUE
#'     ),
#'     list(
#'       formula = Surv(U, D) ~ X,
#'       method = "breslow",
#'       x = TRUE,
#'       y = TRUE
#'     )
#'   ),
#'   predict_fun_list = list(prob_predict.coxph, prob_predict.coxph),
#'   data = dd,
#'   times = seq(1, 5, 0.1),
#'   values = list(X = c(0, 1)),
#'   B = 100,
#'   reference = 0,
#'   contrasts = "difference"
#' )
#' print(x)
#' @export standardize_level
standardize_level <- function(
    fitter_list,
    arguments,
    predict_fun_list,
    data,
    values,
    B = NULL,
    ci_level = 0.95,
    contrasts = NULL,
    reference = NULL,
    seed = NULL,
    times = NULL,
    transforms = NULL,
    progressbar = TRUE) {
  ## Preparation and various checks
  n <- nrow(data)

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
  exposure_names <- colnames(valuesout)

  if (length(fitter_list) != length(predict_fun_list) && length(predict_fun_list) != nrow(valuesout)) {
    stop("need the number fitters, prediction functions and the number of values to be the same")
  }

  fit_outcome <- list()
  for (i in seq_len(length(fitter_list))) {
    fit_outcome[[i]] <- fit_helper(arguments[[i]], fitter_list[[i]], data)
  }
  estimate_fun <- function(valuesout, times, data, fit_outcome, exposure_names) {
    if (is.null(times)) {
      estimates <- rep(NA, nrow(valuesout))
    } else {
      estimates <- matrix(NA, ncol = length(times), nrow = nrow(valuesout))
      colnames(estimates) <- paste("t =", times)
    }

    for (i in seq_len(nrow(valuesout))) {
      # Initialize an empty condition
      subset_condition <- rep(TRUE, nrow(data))

      # Loop through exposure_names and covariate_values
      for (j in seq_along(exposure_names)) {
        var_name <- exposure_names[j]
        var_value <- valuesout[i, j]

        # Update the condition based on the current covariate
        subset_condition <- subset_condition & (data[[var_name]] == var_value)
      }
      data_x <- data[subset_condition, ]
      ## Save the predictions for data_x
      if (is.null(times)) {
        estimates[i] <- mean(predict_fun_list[[i]](object = fit_outcome[[i]], newdata = data_x))
      } else {
        estimates[i, ] <- colMeans(predict_fun_list[[i]](object = fit_outcome[[i]], newdata = data_x, times = times))
      }
    }
    estimates
  }
  estimates <- estimate_fun(valuesout, times, data, fit_outcome, exposure_names)
  estimates_boot <- list()
  if (!is.null(B)) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    pb <- utils::txtProgressBar(
      min = 1,
      max = B,
      style = 3,
      width = 50
    )

    #cat("Bootstrapping... This may take some time... \n")
    for (b in seq_len(B)) {
      if(progressbar) utils::setTxtProgressBar(pb, b)
      data_boot <- data[sample(seq_len(n), replace = TRUE), ]
      fit_outcome_boot <- list()
      for (i in seq_len(length(fitter_list))) {
        fit_outcome_boot[[i]] <- fit_helper(arguments[[i]], fitter_list[[i]], data_boot)
      }
      estimates_boot[[b]] <- estimate_fun(valuesout, times, data_boot, fit_outcome_boot, exposure_names)
    }
  }
  valuesout <- cbind(valuesout, estimates)
  res <- structure(
    list(
      B = B,
      estimates = valuesout,
      fit_outcome = fit_outcome,
      estimates_boot = estimates_boot,
      exposure_names = exposure_names,
      times = times
    )
  )
  format_result_standardize(
    res,
    contrasts,
    reference,
    transforms,
    "plain",
    ci_level,
    "std_custom",
    "summary_standardize"
  )
}

summary_standardize <- function(object, ci_level = 0.95,
                                transform = NULL, contrast = NULL, reference = NULL, ...) {

  B <- object[["B"]]
  est_old_table <- object[["estimates"]]
  n_x_levs <- nrow(est_old_table)
  times <- object[["times"]]
  if (length(times) >= 1) {
    est <- est_old_table[, which(!(colnames(est_old_table) %in% object[["exposure_names"]]))]
  } else {
    est <- est_old_table[["estimates"]]
  }
  if (!is.null(B)) {
    if (length(times) > 1) {
      estimates_boot <- array(NA, c(B, n_x_levs, length(times)))
      for (t in seq_len(length(times))) {
        for (b in seq_len(B)) {
          estimates_boot[b, , t] <- object[["estimates_boot"]][[b]][, t]
        }
      }
    } else {
      estimates_boot <- matrix(unlist(object[["estimates_boot"]]), nrow = B, ncol = n_x_levs, byrow = TRUE)
    }
  }
  if (!is.null(transform)) {
    try_eval_transform <- function(transform, x) {
      est <- tryCatch(
        {
          get(as.character(transform))(x)
        },
        error = function(cond) {
          cond
        }
      )
      if (inherits(est, "simpleError")) {
        stop("transformation failed with error: ", est[["message"]])
      }
      if (any(is.na(est) | is.infinite(est))) {
        stop("transformation failed. Function evaluated at values where it is possibly not well-defined")
      }
      est
    }
    est <- try_eval_transform(transform, est)
    if (!is.null(B)) {
      estimates_boot <- try_eval_transform(transform, estimates_boot)
    }
  }
  if (!is.null(contrast)) {
    if (is.null(reference)) {
      stop("When specifying contrast, reference must be specified as well")
    }
    if (length(object[["exposure_names"]]) > 1L) {
      referencepos <- which(apply(sapply(1:length(object[["exposure_names"]]),
                                         \(i) {
                                           est_old_table[, object[["exposure_names"]]][, i] == reference[i]
                                         }), MARGIN = 1, FUN = \(x) all(x)))
    } else {
      referencepos <- match(reference, est_old_table[, object[["exposure_names"]]])
    }

    if (is.na(referencepos)) {
      stop("reference must be a value in x")
    }
    if (contrast == "difference") {
      contrast_fun <- "-"
    } else if (contrast == "ratio") {
      contrast_fun <- "/"
    } else {
      stop("contrast not supported.")
    }
    if (length(times) > 1L) {
      est <- sweep(est, 2, t(est[referencepos, ]), contrast_fun)
      if (!is.null(B)) {
        for (t_ind in seq_len(length(times))) {
          estimates_boot[, , t_ind] <- get(contrast_fun)(estimates_boot[, , t_ind], estimates_boot[, referencepos, t_ind])
        }
      }
    } else {
      est <- get(contrast_fun)(est, est[referencepos])
      if (!is.null(B)) {
        estimates_boot <- get(contrast_fun)(estimates_boot, estimates_boot[, referencepos])
      }
    }
  }
  if (!is.null(B)) {
    alpha <- 1 - ci_level
    if (length(times) > 1) {
      ci <- apply(estimates_boot, c(2, 3), function(x) c(stats::quantile(x, alpha / 2), stats::quantile(x, 1 - alpha / 2)))
    } else {
      ci <- t(apply(estimates_boot, 2, function(x) c(stats::quantile(x, alpha / 2), stats::quantile(x, 1 - alpha / 2))))
    }
  }

  if (is.factor(reference)) {
    reference <- as.character(reference)
  }
  if (!is.null(contrast)) {
    exposure_name_table <- "Exposure"
    exposure_table <- est_old_table[, object[["exposure_names"]]]
  } else {
    exposure_name_table <- object[["exposure_names"]]
    exposure_table <- est_old_table[, object[["exposure_names"]]]
  }
  if (length(times) > 1) {
    if (is.null(B)) {
      est_table <- cbind(exposure_table, est)
      colnames(est_table) <- c(exposure_name_table, paste0("estimate", " (t=", times, ")"))
    } else {
      res_table <- list()
      for (t_ind in seq_len(length(times))) {
        temp <- cbind(est[, t_ind], t(ci[, , t_ind]))
        colnames(temp) <- c("estimate", "lower", "upper")
        res_table[[t_ind]] <- data.frame(temp)
      }
      est_table <- cbind(exposure_table, do.call("cbind", res_table))
      colnames(est_table)[seq_len(length(exposure_name_table))] <- exposure_name_table
    }
  } else {
    est_table <- data.frame(exposure_table, as.matrix(est, nrow = length(est), ncol = 1L))
    colnames(est_table) <- c(exposure_name_table, "Estimate")
    if (!is.null(B)) {
      ci_boot_df <- data.frame(ci)
      colnames(ci_boot_df) <- paste0(c("lower.", "upper."), ci_level)
      est_table <- cbind(est_table, ci_boot_df)
    }
  }
  out <- c(object, list(
    est_table = est_table, transform = transform,
    contrast = contrast, reference = reference,
    ci_level = ci_level
  ))
  return(out)
}

#' @rdname print
#' @export print.std_custom
#' @export
print.std_custom <- function(x, ...) {
  B <- x[["res"]][["B"]]
  if (!is.null(B)) {
    cat("Number of bootstraps: ", B, "\n")
    cat("Confidence intervals are based on percentile bootstrap confidence intervals \n\n")
  }
  cat("Exposure: ", toString(x[["res"]][["exposure_names"]]), "\n")
  cat("Tables: \n \n")
  for (l in seq_len(length(x[["res_contrast"]]))) {
    temp <- x[["res_contrast"]][[l]]
    if (!is.null(temp[["transform"]])) {
      cat("Transform: ", levels(temp[["transform"]])[[temp[["transform"]]]], "\n")
    }
    if (!is.null(temp[["contrast"]])) {
      cat("Reference level: ", temp[["input"]][["X"]], "=", temp[["reference"]], "\n")
      cat("Contrast: ", temp[["contrast"]], "\n")
    }
    if (is.null(temp[["times"]])) {
      print(temp[["est_table"]], digits = 3L)
    } else {
      if (!is.null(temp[["contrast"]])) {
        len_exposure <- 1
      } else {
        len_exposure <- length(temp[["exposure_names"]])
      }
      for (ti in seq_len(length(temp[["times"]]))) {
        cat("Time: ", temp[["times"]][ti], "\n")
        if (is.null(temp[["B"]])) {
          print(temp[["est_table"]][, c(seq_len(len_exposure), len_exposure + ti)])
        } else {
          print(temp[["est_table"]][, c(seq_len(len_exposure), len_exposure + 3 * ti - c(2, 1, 0))])
        }
        cat("\n")
      }
    }
    cat("\n")
  }
  invisible(x)
}

fit_helper <- function(args, fitter, data) {
  ## try fitting a glm model
  args[["data"]] <- data
  fit <- tryCatch(
    {
      do.call(fitter, args)
    },
    error = function(cond) {
      return(cond)
    }
  )
  if (inherits(fit, "simpleError")) {
    stop("fitter failed with error: ", fit[["message"]])
  } else {
    fit
  }
}
