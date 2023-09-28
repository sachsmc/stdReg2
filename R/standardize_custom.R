#' @title Get standardized estimates using the g-formula with a custom model
#' @inherit standardize_glm
#' @param arguments
#' The arguments to be used in the fitter function as a \code{list}.
#' @param fitter The function to call to fit the data.
#' @param predict_fun The function used to predict the means/probabilities
#' for a new data set on the response level. For survival data,
#' this should be a matrix where each column is the time, and each
#' row the data.
#' @param times For use with survival data. Set to \code{NULL} otherwise.
#' @param B Number of nonparametric bootstrap resamples. Default is \code{NULL} (no bootstrap).
#' @param seed The seed to use with the nonparametric bootstrap.
#' @returns
#' An object of class \code{std_custom}.
#' This is basically a list with components estimates and fit for the outcome model.
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
#'   arguments = list(
#'     formula = Y ~ X * Z,
#'     family = "binomial"
#'   ),
#'   predict_fun = prob_predict.glm,
#'   fitter = "glm",
#'   data = dd,
#'   values = list(X = seq(-1, 1, 0.1)),
#'   B = 100,
#'   references = 0,
#'   contrasts = "difference"
#' )
#' x
#' plot(x)
#' plot(x, reference = 0, contrast = "difference")
#' plot(x, reference = 0, contrast = "difference", plot_ci = FALSE)
#'
#' require(survival)
#' prob_predict.coxph <- function(object, newdata, times){
#'   fit.detail <- suppressWarnings(basehaz(object))
#'   cum.haz <- fit.detail$hazard[sapply(times, function(x) max(which(fit.detail$time <= x)))]
#'   predX <- predict(object = object, newdata = newdata, type = "risk")
#'   res <- matrix(NA, ncol = length(times), nrow = length(predX))
#'   for (ti in seq_len(length(times))){
#'     res[, ti] <- exp(-predX*cum.haz[ti])
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
#'   arguments = list(
#'     formula = Surv(U, D) ~ X + Z + X * Z,
#'     method = "breslow",
#'     x = TRUE,
#'     y = TRUE
#'   ),
#'   fitter = "coxph",
#'   data = dd,
#'   times = 1:5,
#'   predict_fun = prob_predict.coxph,
#'   values = list(X = c(-1, 0, 1)),
#'   B = 100,
#'   references = 0,
#'   contrasts = "difference"
#' )
#' x
#' plot(x)
#' plot(x, reference = 0, contrast = "difference")
#' plot(x, reference = 0, contrast = "difference", plot_ci = FALSE)
#' @export standardize
standardize <- function(arguments,
                        data,
                        fitter,
                        predict_fun,
                        values,
                        B = NULL,
                        ci_level = 0.95,
                        contrasts = NULL,
                        references = NULL,
                        seed = NULL,
                        times = NULL,
                        transforms = NULL) {
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
  exposure <- data[, exposure_names]

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

    cat("Bootstrapping... This may take some time... \n")
    for (b in seq_len(B)) {
      utils::setTxtProgressBar(pb, b)
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
  format_result_standardize(res,
                            contrasts,
                            references,
                            transforms,
                            "plain",
                            ci_level,
                            "std_custom",
                            "summary_standardize")
}

summary_standardize <- function(object, ci_level = 0.95,
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
    reference <- gsub(" ", "", reference, fixed = TRUE)
    if (length(object[["exposure_names"]]) > 1L) {
      est_old_table[["Exposure"]] <- do.call(paste, c(est_old_table[, object[["exposure_names"]]], sep = ","))
    } else {
      est_old_table[["Exposure"]] <- est_old_table[[object[["exposure_names"]]]]
    }

    referencepos <- match(reference, est_old_table[["Exposure"]])
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
    exposure_table <- est_old_table[["Exposure"]]
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
      colnames(est_table)[1:length(exposure_name_table)] <- exposure_name_table
    }
  } else {
    est_table <- data.frame(exposure_table, as.matrix(est, nrow = length(est), ncol = 1L))
    colnames(est_table) <- c(exposure_name_table, "Estimate")
    if (!is.null(B)) {
      ci_boot_df <- data.frame(ci)
      colnames(ci_boot_df) <- paste(c("lower", "upper"), ci_level)
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
    temp <- x[["res_contrast"]][[paste0("V", l)]]
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

summary.plot_help <- function(object,
                              ci_level,
                              transform,
                              contrast,
                              reference,
                              ...) {
  times <- object$input$t
  reference <- ifelse(is.null(reference), "NULL", as.character(reference))
  contrast <- ifelse(is.null(contrast), "NULL", as.character(contrast))
  transform <- ifelse(is.null(transform), "NULL", as.character(transform))
  level_exists <- FALSE
  for (res_cont in object[["res_contrasts"]]) {
    res_cont[["reference"]] <- ifelse(is.null(res_cont[["reference"]]), "NULL", as.character(res_cont[["reference"]]))
    res_cont[["contrast"]] <- ifelse(is.null(res_cont[["contrast"]]), "NULL", as.character(res_cont[["contrast"]]))
    res_cont[["transform"]] <- ifelse(is.null(res_cont[["transform"]]), "NULL", as.character(res_cont[["transform"]]))

    if ((res_cont[["ci_level"]] == ci_level || res_cont[["ci_level"]] == ci_level) &&
      (identical(res_cont[["transform"]], transform)) &&
      (identical(res_cont[["contrast"]], contrast)) &&
      (identical(res_cont[["reference"]], reference))) {
      level_exists <- TRUE
      xlab <- res_cont[["exposure_names"]]
      x <- res_cont[["est_table"]][, 1]
      break
    }
  }
  if (!level_exists) {
    stop("The reference, transform or contrast could not be located in the fitted object.")
  }
  res_table <- res_cont[["est_table"]]
  if (!is.null(times)) {
    res <- list()
    B <- object[["input"]][["B"]]
    for (t_ind in seq_len(length(times))) {
      if (!is.null(B)) {
        res[[t_ind]] <- data.frame(
          Estimate = res_table[, 3 * t_ind - 1],
          Std.Er = NA,
          lower = res_table[, 3 * t_ind],
          upper = res_table[, 3 * t_ind + 1]
        )
      } else {
        res[[t_ind]] <- data.frame(
          Estimate = res_table[, t_ind + 1],
          Std.Er = NA,
          lower = NA,
          upper = NA
        )
      }
    }
    list(est.table = res)
  } else {
    B <- object[["B"]]
    if (is.null(B)) {
      list(est_table = data.frame(res_table[, 1:2], Std.Er. = NA, lower = NA, upper = NA))
    } else {
      list(est_table = data.frame(res_table[, 1:2], Std.Er. = NA, res_table[, 3:4]))
    }
  }
}

#' @rdname plot
#' @export plot.std_custom
#' @export
plot.std_custom <- function(x,
                            plot_ci = TRUE,
                            ci_level = 0.95,
                            transform = NULL,
                            contrast = NULL,
                            reference = NULL, ...) {
  dots <- list(...)
  times <- x[["res"]][["times"]]
  B <- x[["res"]][["B"]]

  if (length(x[["res"]][["exposure_names"]]) > 1L) {
    stop("cannot do plot with multiple exposures")
  }

  if (!is.null(times)) {
    obj <- list(res = list(
      input = list(
        valuesout = x[["res"]][["estimates"]][1],
        times = times,
        B = B
      ),
      res_contrasts = x[["res_contrast"]])
    )
    class(obj) <- "plot_help"
    plot.std_surv(
      x = obj,
      plot_ci = !is.null(B) && plot_ci,
      ci_level = ci_level,
      transform = transform,
      reference = reference,
      contrast = contrast,
      summary_fun = "summary.plot_help"
    )
  } else {
    temp <- x
    temp[["res"]][["res_contrasts"]] <- x[["res_contrast"]]
    class(temp[["res"]]) <- "plot_help"
    ## needs summary function to be rewriten
    plot.std_glm(
      x = temp,
      plot_ci = !is.null(B) && plot_ci,
      ci_level = ci_level,
      transform = transform,
      reference = reference,
      contrast = contrast,
      summary_fun = "summary.plot_help"
    )
  }
}
