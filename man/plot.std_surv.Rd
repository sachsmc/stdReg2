% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxph_methods.R
\name{plot.std_surv}
\alias{plot.std_surv}
\title{Plots regression standardization fit}
\usage{
\method{plot}{std_surv}(
  x,
  plot_ci = TRUE,
  ci_type = "plain",
  ci_level = 0.95,
  transform = NULL,
  contrast = NULL,
  reference = NULL,
  legendpos = "bottomleft",
  summary_fun = "summary_std_coxph",
  ...
)
}
\arguments{
\item{x}{An object of class \code{"std_surv"}.}

\item{plot_ci}{if \code{TRUE}, add the confidence intervals to the plot.}

\item{ci_type}{A string, indicating the type of confidence intervals. Either "plain", which
gives untransformed intervals, or "log", which gives log-transformed intervals.}

\item{ci_level}{Coverage probability of confidence intervals.}

\item{transform}{If set to \code{"log"}, \code{"logit"}, or \code{"odds"}, the standardized
mean \eqn{\theta(x)} is transformed into \eqn{\psi(x)=\log\{\theta(x)\}},
\eqn{\psi(x)=\log[\theta(x)/\{1-\theta(x)\}]}, or
\eqn{\psi(x)=\theta(x)/\{1-\theta(x)\}}, respectively. If left unspecified,
\eqn{\psi(x)=\theta(x)}.}

\item{contrast}{If set to \code{"difference"} or \code{"ratio"}, then \eqn{\psi(x)-\psi(x_0)}
or \eqn{\psi(x) / \psi(x_0)} are constructed, where \eqn{x_0} is a reference
level specified by the \code{reference} argument.
If not \code{NULL}, a doubly robust estimator of the standardized estimator is used.}

\item{reference}{If \code{contrast} is specified, the desired reference level.}

\item{legendpos}{position of the legend; see \link[graphics]{legend}.}

\item{summary_fun}{For internal use only. Do not change.}

\item{\dots}{Unused.}
}
\value{
None. Creates a plot as a side effect
}
\description{
This is a \code{plot} method for class \code{"std_surv"}.
}
