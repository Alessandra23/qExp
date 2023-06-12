#' \eqn{g(\theta)} function
#'
#' Function to calculate \eqn{g(\theta)}.
#'
#' @param u Tuning parameter
#' @param v Tuning parameter
#' @param theta Parameter \eqn{\theta} of qExponential function
#' @export


g.theta <- function(theta, v, u = -v) {
  (theta + 1)^(u - v) * beta(v + 1, theta - v + 1)^u * beta(u + 1, theta - u + 1)^(-v)
}

#' \eqn{g(\theta)} function inverse
#'
#' Calculate the inverse of the function \eqn{g(\theta)}.
#'
#' @param u Tuning parameter
#' @param v Tuning parameter
#' @param theta Parameter \eqn{\theta} of qExponential function
#'
#' @import  GoFKernel
#' @export
g.theta.inv <- function(value, lower = 0, upper = 10000, v, u = -v) {
  g <- function(theta) {
    (theta + 1)^(u - v) * beta(v + 1, theta - v + 1)^u * beta(u + 1, theta - u + 1)^(-v)
  }
  g.inverse <- GoFKernel::inverse(g, lower = lower, upper = upper)
  inverse.theta <- g.inverse(value)
  return(inverse.theta)
}

#' Limits of \eqn{g(\theta)} function
#'
#' @param u Tuning parameter
#' @param v Tuning parameter
#' @return A list of two elements: limit of \eqn{g(\theta)} when \eqn{\theta} goes to 0 and limit of \eqn{g(\theta)} when \eqn{\infty}
#' @export
limits.g <- function(v,u){
  lim.inf <- ((gamma(v) * v)^u / ((gamma(u) * u)^v)) * ((-gamma(-v) * v)^u / ((-gamma(-u) * u)^v))
  # lim.sup <- 2^(-(u+v)/2) * exp(((log(2)+log(pi))*(u+v))/2)*pi^(-(u+v)/2) * (gamma(v+1)^u/gamma(u+1)^v)
  lim.sup <-  (gamma(v+1)^u/gamma(u+1)^v)
  return(list(lim.inf = lim.inf,
              lim.sup = lim.sup))
}


#'
#'
#'
#' @export

f.vu <- function(par, mu = 1, theta = 1/9, prob = 0.97){
  v <- par[1]
  u <- par[2]

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
  lim.sup <- (gamma(v+1)^u/gamma(u+1)^v)

  return(ceiling((qnorm(prob)*sqrt(gamma2)/(lim.sup-Esp.Tn))^2))
}

