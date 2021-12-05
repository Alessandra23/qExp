#' Theoretical quantities MFMM
#'
#' Theoretical quantities for the mfmm
#' @param mu true value of mu
#' @param theta true value of thete
#' @param v Tuning parameter
#' @param u Tuning parameter
#' @export
mfmm.theo <- function(mu, theta, u, v) {
  ev <- (mu * theta)^v * (theta + 1) * beta(v + 1, theta + 1 - v)
  eu <- (mu * theta)^u * (theta + 1) * beta(u + 1, theta + 1 - u)
  sig.v <- (mu * theta)^(2 * v) * (theta + 1) * (beta(2 * v + 1, theta + 1 - 2 * v) - (theta + 1) * (beta(v + 1, theta + 1 - v)^2))
  sig.u <- (mu * theta)^(2 * u) * (theta + 1) * (beta(2 * u + 1, theta + 1 - 2 * u) - (theta + 1) * (beta(u + 1, theta + 1 - u)^2))
  sig.uv <- ((mu * theta)^(u + v)) * (theta + 1) * (beta(u + v + 1, theta + 1 - u - v) - (theta + 1) * beta(u + 1, theta + 1 - u) * beta(v + 1, theta + 1 - v))
  p1 <- (u^2) * sig.v * ev^(2 * u - 2) * eu^(-2 * v)
  p2 <- u * v * ev^(2 * u - 1) * eu^(-2 * v - 1) * sig.uv
  p3 <- (v^2) * sig.u * ev^(2 * u) * eu^(-2 * v - 2)
  gamma2 <- p1 - 2 * p2 + p3
  Esp.Tn <- (ev^u) / (eu^v)

  return(list(
    Esp.Tn = Esp.Tn,
    gamma2 = gamma2
  ))
}


#' Observed quantities MFMM
#'
#' Observed quantities for the MFMM
#'
#' @import numDeriv
#' @param n sample size
#' @param mu true value of mu
#' @param theta true value of theta
#' @param v Tuning parameter
#' @param u Tuning parameter
#' @param samples matrix os samples
#' @export
mfmm.samples <- function(n, mu, theta, v, u = -v, samples, d = 1e-10, lower = 0, upper = 10000) {
  theo <- mfmm.theo(mu, theta, u, v)
  Esp.Tn <- theo$Esp.Tn
  gamma2 <- theo$gamma2

  g <- function(theta) {
    g.theta(theta, v = v, u = u)
  }
  g.inverse <- inverse(g, lower = lower, upper = upper)

  inv.esp <- g.inverse(Esp.Tn)
  der.inv <- numDeriv::grad(func = g.inverse, x = Esp.Tn, method.args = list(eps = 1e-12, d = d, r = 6))
  # pracma::grad(g.inverse, Esp.Tn)

  lim.inf <- ((gamma(v) * v)^u / ((gamma(u) * u)^v)) * ((-gamma(-v) * v)^u / ((-gamma(-u) * u)^v))
  lim.sup <- (10000 + 1)^(u - v) * beta(v + 1, 10000 - v + 1)^u * beta(u + 1, 10000 - u + 1)^(-v)

  mu.uv <- (rowMeans(samples^v)^u) / (rowMeans(samples^u)^v)

  # Select sample
  ifelse(v > 0, samp.cond <- mu.uv[mu.uv > lim.inf & mu.uv < lim.sup],
    samp.cond <- mu.uv[mu.uv < lim.inf & mu.uv > lim.sup]
  )

  # estimated theta based on selected sample
  theta.hat <- sapply(samp.cond, g.theta.inv, v = v, u = u)
  # variance obtained by the delta method for the inverse of the function
  kappa2 <- gamma2 * (der.inv)^2
  # standardizing the estimated theta
  theta.hat.pad <- (sqrt(n) * (theta.hat - theta)) / sqrt(kappa2)
  # q estimates
  q.hat <- (3 + theta.hat) / (2 + theta.hat)
  # Asymptotic variance of q
  var.qhat <- (inv.esp + 2)^(-4) * kappa2
  # true value of q
  q.real <- (3 + theta) / (2 + theta)
  # standardizing q estimated
  q.hat.pad <- (sqrt(n) * (q.hat - q.real)) / sqrt(var.qhat)

  # d.q.theta <- data.frame(theta.hat,theta.hat.pad,q.hat,q.hat.pad)

  return(list(
    theta.hat = theta.hat, theta.hat.pad = theta.hat.pad,
    q.hat = q.hat, q.hat.pad = q.hat.pad,
    kappa2 = kappa2, lim.sup = lim.sup
  ))
}
