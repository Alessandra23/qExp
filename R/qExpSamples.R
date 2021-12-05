#' Samples generator
#'
#' Generates samples in the parameterization in terms of \eqn{\mu} and \eqn{\theta} (Patriota, 2012)
#'
#'
#' @param N Number of samples
#' @param n Sample size of each sample
#' @param theta True value of theta
#' @param mu True value of mu
#'
#' @return A matrix of order N x n, where each row is a sample of size n.
#'
#' @export
#'
qexp.samples <- function(N,n,theta,mu){
  samples <- t(replicate(N, rtsal(n,shape = theta+1, scale = mu*theta)))
  return(samples)
}

