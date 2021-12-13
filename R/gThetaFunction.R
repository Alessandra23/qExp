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
#' @importFrom  GoFKernel inverse
#' @export
g.theta.inv <- function(value, lower = 0, upper = 10000, v, u = -v) {
  g <- function(theta) {
    (theta + 1)^(u - v) * beta(v + 1, theta - v + 1)^u * beta(u + 1, theta - u + 1)^(-v)
  }
  g.inverse <- inverse(g, lower = lower, upper = upper)
  inverse.theta <- g.inverse(g(value))
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
  lim.sup <- 2^(-(u+v)/2) * exp(((log(2)+log(pi))*(u+v))/2)*pi^(-(u+v)/2) * (gamma(v+1)^u/gamma(u+1)^v)
  return(list(lim.inf = lim.inf,
              lim.sup = lim.sup))
}
