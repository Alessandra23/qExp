#' q-Exponential generator

#' Distribution-related (per R standards):
#' dtsal				Probability density
#' ptsal				Cumulative probability
#' qtsal				Quantiles
#' rtsal				Random variate generation
#' Parameter Conversion:
#' tsal.shape.from.q		Get shape parameter from q
#' tsal.scale.from.qk		Get scale parameter from q, kappa
#' tsal.q.from.shape		Get q from from shape parameter
#' tsal.kappa.from.ss		Get kappa from shape, scale
#' tsal.ss.from.qk		Get shape, scale from q, kappa (as pairs)
#' tsal.qk.from.ss		Get q, kappa from shape, scale (as pairs)


#################################################################
########        PARAMETER CONVERSION FUNCTIONS          #########
#################################################################

#' Users will generally have little reason to invoke these functions
#'
#' Calculate generalized-Pareto shape parameter from q parameter
#' @param q a real value
#' @return a real value
#' @export
tsal.shape.from.q <- function(q) {
  shape <- -1 / (1 - q)
  return(shape)
}

#' Calculate generalized-Pareto scale parameter from q and kappa parameters
#' Input: two real values
#' Output: a real value
#' @export
tsal.scale.from.qk <- function(q, kappa) {
  shape <- tsal.shape.from.q(q)
  scale <- shape * kappa
  return(scale)
}

#' Calculate q parameter from shape parameter
#' Input: a real value
#' Output: a real value
#' @export
tsal.q.from.shape <- function(shape) {
  q <- 1 + 1 / shape
  return(q)
}

#' Calculate kappa parameter from shape and scale parameters
#' Input: two real values
#' Output: a real value
#' @export
tsal.kappa.from.ss <- function(shape, scale) {
  kappa <- scale / shape
  return(kappa)
}

#' Calculate shape & scale parameters from q & kappa parameters
#' Input: two real values
#' Output: vector of two real values
#' @export
tsal.ss.from.qk <- function(q, kappa) {
  ss <- c(tsal.shape.from.q(q), tsal.scale.from.qk(q, kappa))
  return(ss)
}

#' Calculate q & kappa parameters from shape & scale parameters
#' Input: two real values
#' Output: vector of two real values
#' @export
tsal.qk.from.ss <- function(shape, scale) {
  qk <- c(tsal.q.from.shape(shape), tsal.kappa.from.ss(shape, scale))
  return(qk)
}

#################################################################
########       DISTRIBUTION-RELATED FUNCTIONS           #########
#################################################################


#' Calculate the probability density
#' Input: vector of data values, distributional parameters, left-censoring
#'        threshold, log flag
#' Output: vector of (log) densities
#' @export
dtsal <- function(x, shape = 1, scale = 1, q = tsal.q.from.shape(shape),
                  kappa = tsal.kappa.from.ss(shape, scale), xmin = 0,
                  log = FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q, kappa)
  shape <- ss[1]
  scale <- ss[2]
  # Under censoring, pass off to the tail version
  if (xmin > 0) {
    return(dtsal.tail(x, shape, scale, q, kappa, xmin, log))
  }
  z <- 1 + x / scale
  if (log) {
    d <- -(shape + 1) * log(z) + log(shape / scale)
  } else {
    d <- (shape / scale) * (z^(-shape - 1))
  }
  return(d)
}

#' Calculate the cumulative distribution function
#' Input: vector of data values, distributional parameters, left-censoring
#'        threshold, usual flags
#' Output: vector of (log) probabilities
#' @export
ptsal <- function(x, shape = 1, scale = 1, q = tsal.q.from.shape(shape),
                  kappa = tsal.kappa.from.ss(shape, scale), xmin = 0,
                  lower.tail = TRUE, log.p = FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q, kappa)
  shape <- ss[1]
  scale <- ss[2]
  if (xmin > 0) {
    return(ptsal.tail(
      x, shape, scale, q, kappa, xmin, lower.tail,
      log.p
    ))
  }
  z <- 1 + x / scale
  if ((log.p) && (!lower.tail)) {
    p <- -shape * log(z)
  }
  if ((log.p) && (lower.tail)) {
    p <- log(1 - (z^(-shape)))
  }
  if ((!log.p) && (!lower.tail)) {
    p <- z^(-shape)
  }
  if ((!log.p) && (lower.tail)) {
    p <- 1 - z^(-shape)
  }
  return(p)
}

#' Calculate quantiles
#' Input: vector of p-values, distributional parameters, left-censoring
#'        threshold, usual flags
#' Output: vector of quantile locations
#' @export
qtsal <- function(p, shape = 1, scale = 1, q = tsal.q.from.shape(shape),
                  kappa = tsal.kappa.from.ss(shape, scale), xmin = 0,
                  lower.tail = TRUE, log.p = FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q, kappa)
  shape <- ss[1]
  scale <- ss[2]
  if (xmin > 0) {
    return(qtsal.tail(
      p, shape, scale, q, kappa, xmin, lower.tail,
      log.p
    ))
  }
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail) {
    p <- 1 - p
  }
  # The upper quantile function is given by
  # (scale)(p^{-1/shape} - 1) = x
  quantiles <- scale * (-1 + (p^(-1 / shape)))
  return(quantiles)
}

#' Generate random variates
#' Input: integer length, distributional parameters, left-censoring threshold
#' Output: vector of reals
#' @export
rtsal <- function(n, shape = 1, scale = 1, q = tsal.q.from.shape(shape),
                  kappa = tsal.kappa.from.ss(shape, scale), xmin = 0) {
  # If we have both shape/scale and q/kappa parameters, the latter over-ride.
  ss <- tsal.ss.from.qk(q, kappa)
  shape <- ss[1]
  scale <- ss[2]
  if (xmin > 0) {
    return(rtsal.tail(n, shape, scale, q, kappa, xmin))
  }
  # Apply the transformation method
  ru <- runif(n)
  r <- qtsal(ru, shape, scale)
  return(r)
}
