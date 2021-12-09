#' Log of likelihood qexp
#' @param param vector of paramters \eqn{(\mu,\theta)}
#' @param n sample size
#' @export
log.qexp <- function(param, n) {
  mu <- param[1]
  theta <- param[2]
  return(-(n * log((theta + 1) / (theta * mu)) - (theta + 2) * sum(log(1 + x / (mu * theta)))))
}

#' MLE estimation
#' @param sample Sample of size n
#' @export
est.mle <- function(sample) {
  mu.est <- mean(sample)
  theta.est <- mean(sample^0.1) / (mean(sample)^0.1)
  est.mu.theta <- optim(par = c(mu.est, theta.est), log.qexp, n = length(sample))
  return(est.mu.theta$par)
}

#' Expected Fisher Information
#' @param mu True value of mu
#' @param theta True value of theta
#' @export
K.mu.theta <- function(mu, theta) {
  KF <- matrix(0, 2, 2)
  KF[1, 1] <- ((theta + 2)^2) / (mu^2 * (theta^2 + 5 * theta + 6))
  KF[2, 2] <- (theta^2 + theta + 2) / (theta^2 * (theta + 1)^2 * (theta^2 + 5 * theta + 6))
  KF[1, 2] <- 2 / (mu * theta * (theta^2 + 5 * theta + 6))
  KF[2, 1] <- KF[1, 2]
  return(KF)
}

#' Theoretical variances
#' @param mu True value of mu
#' @param theta True value of theta
#' @export
var.teo.mle <- function(mu, theta) {
  var.mu <- solve(K.mu.theta(c(mu, theta)))[1, 1]
  var.theta <- solve(K.mu.theta(c(mu, theta)))[2, 2]
  return(list(var.mu = var.mu, var.theta = var.theta))
}






#' Sample quantities MLE
#'
#' Theoretical quantities for the mle
#' @param samples matrix of samples Nxn
#' @param mu true value of mu
#' @param theta true value of theta
#' @export
mle.samples <- function(samples, mu, theta) {
  n <- ncol(samples)
  est.mle <- t(apply(samples, 1, est.mle))
  est.mle.mu <- est.mle[, 1]
  est.mle.theta <- est.mle[, 2]
  est.mle.q <- (3 + est.mle.theta) / (2 + est.mle.theta)
  q.real <- (theta + 3) / (theta + 2)
  var.mu <- solve(K.mu.theta(mu, theta))[1, 1] ## Of mu
  var.theta <- solve(K.mu.theta(mu, theta))[2, 2] ## Of theta
  mu.pad.mle <- sqrt(n) * (est.mle.mu - mu) / sqrt(var.mu)
  theta.pad.mle <- sqrt(n) * (est.mle.theta - theta) / sqrt(var.theta)
  q.pad.mle <- sqrt(n) * (est.mle.q - q.real) / sqrt(var.theta * ((theta + 2)^(-4)))

  return(list(
    est.mle.mu = est.mle.mu, est.mle.theta = est.mle.theta,
    est.mle.q = est.mle.q, mu.pad.mle = mu.pad.mle,
    theta.pad.mle = theta.pad.mle, q.pad.mle = q.pad.mle
  ))
}
