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
g.theta.inv <- function(theta, lower = 0, upper = 1000, v, u = -v) {
  g <- function(theta) {
    (theta + 1)^(u - v) * beta(v + 1, theta - v + 1)^u * beta(u + 1, theta - u + 1)^(-v)
  }
  g.inverse <- inverse(g, lower = lower, upper = upper)
  inverse.theta <- g.inverse(g(theta))
  return(inverse.theta)
}
