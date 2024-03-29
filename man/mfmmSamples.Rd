% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mfmmQuantities.R
\name{mfmmSamples}
\alias{mfmmSamples}
\title{Observed quantities MFMM}
\usage{
mfmmSamples(n, mu, theta, v, u, samples, d = 1e-10, lower = 0, upper = 10000)
}
\arguments{
\item{n}{sample size}

\item{mu}{true value of mu}

\item{theta}{true value of theta}

\item{v}{Tuning parameter}

\item{u}{Tuning parameter}

\item{samples}{matrix os samples}

\item{d}{parameter to calculate the gradient}

\item{lower}{lower to calculate the inverse of g}

\item{upper}{upper to calculate the inverse og g}
}
\description{
Observed quantities for the MFMM
}
